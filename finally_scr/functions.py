import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from scipy.sparse import coo_matrix, csr_matrix
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, auc, confusion_matrix, accuracy_score, precision_score, recall_score


def load_data(chrom_nodes_dict, chrom_edges, shuffle_data=True):
    #节点字典直接获取节点ID
    enhancer_ids = chrom_nodes_dict['enhancer']['ids']
    promoter_ids = chrom_nodes_dict['promoter']['ids']
    enhancer_features = chrom_nodes_dict['enhancer']['features']
    promoter_features = chrom_nodes_dict['promoter']['features']

    enhancer_index = {id: i for i, id in enumerate(enhancer_ids)}
    promoter_index = {id: i + len(enhancer_ids) for i, id in enumerate(promoter_ids)}
    edges = []
    for idx, row in chrom_edges.iterrows():
        enhancer_idx = enhancer_index.get(row['ccRE_ID'])
        promoter_idx = promoter_index.get(row['Gene_ID'])
        if enhancer_idx is None or promoter_idx is None:
            print(f"Missing index for ccRE_ID {row['ccRE_ID']} or Gene_ID {row['Gene_ID']}")
        else:
            edges.append([enhancer_idx, promoter_idx])
    # edges = np.array(edges).T   #(2, E)
    edges = np.array(edges)   #(E, 2)
    
    # 加载标签和距离信息
    labels = chrom_edges['label'].values.astype(int)   #[E]
    distances = chrom_edges['distance'].values.astype(float)  #[E]

    # print(edges.shape,labels.shape,distances.shape)
    # 随机打乱数据（如果需要）
    if shuffle_data:
        indices = np.random.permutation(edges.shape[0])  
        edges = edges[indices]
        labels = labels[indices]
        distances = distances[indices]

    return {
        'enhancer_features': enhancer_features,
        'promoter_features': promoter_features,
        'edges': edges,
        'labels': labels,
        'distances': distances
    }


def create_data_objects(chrom_data, is_training=True):
    data_list = []
    enhancer_features = torch.tensor(chrom_data['enhancer_features'], dtype=torch.float)
    promoter_features = torch.tensor(chrom_data['promoter_features'], dtype=torch.float)
    edge_index = torch.tensor(chrom_data['edges'].T, dtype=torch.long)       
    num_nodes = enhancer_features.size(0) + promoter_features.size(0)
    # 计算节点的度
    node_degree = torch.zeros(num_nodes, dtype=torch.long)  
    ones = torch.ones(edge_index.size(1), dtype=torch.long)  
    node_degree = node_degree.scatter_add(0, edge_index[0], ones) 
    node_degree = node_degree.scatter_add(0, edge_index[1], ones) 
    # 转换为float并标准化节点度
    node_degree = node_degree.float().unsqueeze(1)
    node_degree = (node_degree - node_degree.mean()) / (node_degree.std() + 1e-5)  # 添加一个小的常数防止除以0

    # 边属性：距离和标签
    edge_distances = torch.tensor(chrom_data['distances'], dtype=torch.float).unsqueeze(1) #增加维度，[E, 1]
    edge_labels = torch.tensor(chrom_data['labels'], dtype=torch.long).unsqueeze(1)
    
    #将标签作为边特征，但是引入掩码技术，在训练的时候使用，但是测试的时候掩盖
    if is_training:
        edge_attr = torch.cat([edge_distances, edge_labels.float()], dim=1)  #[E, 2]
    else:
        edge_attr = torch.cat([edge_distances, torch.zeros_like(edge_labels.float())], dim=1)

    if torch.isnan(edge_labels).any():
        print("NaN detected in labels")
    if torch.isnan(enhancer_features).any():
        print("NaN detected in enhancer_features")
    if torch.isnan(promoter_features).any():
        print("NaN detected in promoter_features")
    if torch.isnan(edge_attr).any():
        print("NaN detected in edge_attr")
    if torch.isnan(node_degree).any():
        print("NaN detected in 标准化后的node_degree")

    # 创建 Data 对象
    data_object = Data(x_e=enhancer_features, x_p=promoter_features, edge_index=edge_index,
                       edge_attr=edge_attr, y=edge_labels, num_nodes=num_nodes, node_degree=node_degree)
    data_list.append(data_object)
    return data_list


def load_subgraph_loader(nodes_dict, edges_df, geometric_batch_size, is_training=True):
    chrom_dataloaders = {}
    unique_chromosomes = pd.unique(edges_df['chr'])

    for chrom in unique_chromosomes:
        chrom_edges = edges_df[edges_df['chr'] == chrom]

        # 使用numpy布尔索引筛选当前染色体的enhancer和promoter
        enhancer_mask = nodes_dict['enhancer']['chr'] == chrom
        promoter_mask = nodes_dict['promoter']['chr'] == chrom

        # 直接使用筛选后的数组，无需转换为列表
        chrom_nodes_dict = {
            'enhancer': {
                'ids': nodes_dict['enhancer']['ids'][enhancer_mask],
                'features': nodes_dict['enhancer']['features'][enhancer_mask]
            },
            'promoter': {
                'ids': nodes_dict['promoter']['ids'][promoter_mask],
                'features': nodes_dict['promoter']['features'][promoter_mask]
            }
        }

        # 加载当前染色体的数据
        chrom_data = load_data(chrom_nodes_dict, chrom_edges, shuffle_data=is_training)
        data_list = create_data_objects(chrom_data, is_training)
        chrom_dataloaders[chrom] = DataLoader(data_list, batch_size=geometric_batch_size, shuffle=True)

    return chrom_dataloaders



def train_model(model, dataloaders, criterion, optimizer, threshold, device):
    model.train()  # 确保模型处于训练模式
    total_loss = 0.0
    correct_predictions = 0
    total_predictions = 0

    # 遍历每个染色体的 DataLoader
    for chrom, loader in dataloaders.items():
        for data in loader:
            # 将数据发送到设备
            data = data.to(device)   
            optimizer.zero_grad() # 清空过往梯度
            # 模型输出
            outputs = model(data)
            # 计算损失，注意 data.y 已经是 labels
            # loss = criterion(outputs.squeeze(), data.y.float())
            loss = criterion(outputs.squeeze(), data.y.float().squeeze())
            loss.backward()
            optimizer.step()

            total_loss += loss.item() * data.num_graphs
            predictions = (outputs >= threshold).float()
            correct_predictions += torch.sum(predictions == data.y.float()).item()
            total_predictions += data.y.size(0)
    
    avg_loss = total_loss / total_predictions
    acc = correct_predictions / total_predictions
    
    return avg_loss, acc


def evaluate_model(model, dataloader, criterion, threshold, device):
    model.eval()
    y_true, y_pred = [], []
    total_loss = 0.0
    total_samples = 0

    with torch.no_grad():
        for chrom, loader in dataloader.items():
            for data in loader:
                data = data.to(device)
                outputs = model(data)
                # loss = criterion(outputs, data.y.float())
                loss = criterion(outputs.squeeze(), data.y.float().squeeze())
                total_loss += loss.item() * data.num_graphs

                # 保存真实标签和预测概率
                predictions = torch.sigmoid(outputs).squeeze()
                y_pred.extend(predictions.cpu().numpy())
                y_true.extend(data.y.cpu().numpy())

                total_samples += data.num_graphs

    avg_loss = total_loss / total_samples

    # 计算其他性能指标
    auc = roc_auc_score(y_true, y_pred)
    aupr = average_precision_score(y_true, y_pred)
    y_pred_binary = (np.array(y_pred) >= threshold).astype(int)
    acc = accuracy_score(y_true, y_pred_binary)

    return avg_loss, auc, aupr, acc


def hinge_loss(outputs, labels):
    # 将标签从{0, 1}转换为{-1, 1}，以适应铰链损失的标准形式
    labels_transformed = 2 * labels - 1
    
    # 计算铰链损失，即确保正确标签的输出足够远离决策边界,# 铰链损失定义为 max(0, 1 - y_true * y_pred)
    hinge_loss_value = torch.mean(torch.clamp(1 - labels_transformed * outputs, min=0))
    return hinge_loss_value

class CombinedLoss(torch.nn.Module):
    def __init__(self, weight_bce=0.5):
        super(CombinedLoss, self).__init__()
        self.bce_loss = torch.nn.BCEWithLogitsLoss()
        self.weight_bce = weight_bce
        self.weight_hinge = 1 - weight_bce

    def forward(self, outputs, labels):
        # 计算交叉熵损失
        bce_loss = self.bce_loss(outputs, labels)
        # 计算铰链损失
        hinge_loss_value = hinge_loss(outputs, labels)
        # 结合两种损失
        combined_loss = self.weight_bce * bce_loss + self.weight_hinge * hinge_loss_value
        return combined_loss



# def standardize_features(features):
#     mean = np.mean(features, axis=0)
#     std = np.std(features, axis=0)
#     standardized_features = (features - mean) / std
#     return standardized_features


# def load_subgraph_data(node_features_file, edges_file, device):
#     
#     # 加载数据
#     adj, features, edges_dict, labels_dict, distances_dict, chrom = load_data(node_features_file, edges_file)
#     chrom_list = np.unique(chrom)

#     # 创建存储每个子图数据的字典
#     subgraph_edge_index = {}  
#     subgraph_features = {}
#     subgraph_edge_features = {}
#     data_dict = {}

#     # 加载并转换每个子图的数据
#     for chrom in chrom_list:
#         # 存储边索引
#         subgraph_edge_index[chrom] = torch.transpose(torch.tensor(edges_dict[chrom], dtype=torch.long), 0, 1).to(device)
#         # 存储节点特征
#         subgraph_features[chrom] = torch.FloatTensor(features[chrom].toarray()).to(device) if isinstance(features[chrom], csr_matrix) else torch.FloatTensor(features[chrom]).to(device)
#         # 标准化并存储边的特征
#         distances_dict[chrom] = standardize_features(distances_dict[chrom])
#         subgraph_edge_features[chrom] = torch.FloatTensor(distances_dict[chrom]).to(device)
#         # 存储标签
#         labels = torch.FloatTensor(labels_dict[chrom]).to(device)
#         # 组装数据
#         data_dict[chrom] = [subgraph_features[chrom], subgraph_edge_index[chrom], subgraph_edge_features[chrom], labels]

#     node_feature_dim = subgraph_features['chr1'].shape[1]
    
#     return data_dict, chrom_list, node_feature_dim


# def train_model(model, dataloader, criterion, optimizer, threshold, device):
#     model.train()  # 确保模型处于训练模式
#     total_loss = 0.0
#     correct_predictions = 0
#     total_predictions = 0 

#     for data in dataloader:
#         # 将数据发送到设备
#         data = data.to(device)
        
#         optimizer.zero_grad()
#         # 模型输出
#         outputs = model(data)
#         # 计算损失，注意 data.y 已经是 labels
#         loss = criterion(outputs.squeeze(), data.y.float())
#         loss.backward()
#         optimizer.step()

#         total_loss += loss.item() * data.num_graphs
#         predictions = (outputs >= threshold).float()
#         correct_predictions += torch.sum(predictions == data.y.float()).item()
#         total_predictions += data.y.size(0)
    
#     avg_loss = total_loss / total_predictions
#     acc = correct_predictions / total_predictions
    
#     return avg_loss, acc

# def evaluate_model(model, dataloader, criterion, threshold, device):
#    
#     model.eval()
#     y_true, y_pred = [], []
#     total_loss = 0.0
#     total_samples = 0

#     with torch.no_grad():
#         for data in dataloader:
#             data = data.to(device)
#             outputs = model(data)
#             loss = criterion(outputs, data.y.float())
#             total_loss += loss.item() * data.num_graphs

#             # 保存真实标签和预测概率
#             predictions = torch.sigmoid(outputs).squeeze()
#             y_pred.extend(predictions.cpu().numpy())
#             y_true.extend(data.y.cpu().numpy())

#             total_samples += data.num_graphs

#     avg_loss = total_loss / total_samples

#     # 计算其他性能指标
#     auc = roc_auc_score(y_true, y_pred)
#     aupr = average_precision_score(y_true, y_pred)
#     y_pred_binary = (np.array(y_pred) >= threshold).astype(int)
#     acc = accuracy_score(y_true, y_pred_binary)

#     return avg_loss, auc, aupr, acc


# def evaluate_model(model, dataloader, criterion, threshold, device):
#     """
#     评估模型在整个验证集上的性能，并返回平均损失和其他性能指标。

#     参数:
#     - model: 训练好的模型。
#     - dataloader: 验证集的 DataLoader。
#     - criterion: 损失函数。
#     - threshold: 分类阈值，用于将概率转换为二进制标签。
#     - device: 设备，如"cuda"或"cpu"。
#     """
#     model.eval()
#     y_true, y_pred = [], []
#     total_loss = 0.0
#     total_samples = 0

#     with torch.no_grad():
#         for data in dataloader:
#             data = data.to(device)
#             outputs = model(data.x, data.edge_index, data.edge_attr)
#             loss = criterion(outputs, data.y.float())
#             total_loss += loss.item() * data.num_graphs

#             # 保存真实标签和预测概率
#             predictions = torch.sigmoid(outputs).squeeze()
#             y_pred.extend(predictions.cpu().numpy())
#             y_true.extend(data.y.cpu().numpy())

#             total_samples += data.num_graphs

#     avg_loss = total_loss / total_samples

#     # 计算其他性能指标
#     auc = roc_auc_score(y_true, y_pred)
#     aupr = average_precision_score(y_true, y_pred)
#     y_pred_binary = (np.array(y_pred) >= threshold).astype(int)
#     acc = accuracy_score(y_true, y_pred_binary)
#     precision = precision_score(y_true, y_pred_binary, zero_division=0)
#     recall = recall_score(y_true, y_pred_binary, zero_division=0)
#     cm = confusion_matrix(y_true, y_pred_binary)

#     return avg_loss, auc, aupr, acc, precision, recall, cm
