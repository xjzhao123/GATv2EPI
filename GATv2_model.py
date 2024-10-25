import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATv2Conv

class FeatureAugmentation(nn.Module):
    """
    特征增强模块，用于生成增强的节点特征，同时尽量保持协方差的一致性。

    参数:
    num_features (int): 输入特征的维数。
    epsilon (float): 控制噪声强度的超参数，用于确定噪声的方差。
    """
    def __init__(self, num_features, epsilon=0.01):
        super(FeatureAugmentation, self).__init__()
        self.epsilon = epsilon
        # 初始化变换矩阵 P，用于特征变换
        self.P = nn.Parameter(torch.randn((num_features, num_features)))

    def forward(self, x):
        """
        前向传播逻辑，应用特征变换并添加噪声。

        参数:
        x (Tensor): 输入特征矩阵，形状为 [batch_size, num_features]

        返回:
        x_augmented (Tensor): 增强后的特征矩阵。
        """
        # 应用变换矩阵 P 以变换特征
        x_transformed = torch.mm(x, self.P)

        # 添加服从高斯分布的噪声，均值为0，标准差为 epsilon
        noise = torch.randn_like(x_transformed) * self.epsilon

        # 结合变换后的特征和噪声得到增强特征
        x_augmented = x_transformed + noise
        return x_augmented

        # 注意：这里没有直接实施控制协方差一致性的约束
        # 为实现该功能，需要在训练过程中检查和调整 P 和 epsilon

def one_hot_encode(x, type, device):
    """
    生成独热编码向量。

    """
    # 使用独热编码标识节点类型
    if type == 'e':
        one_hot = torch.zeros(x.size(0), 2, device=device)  # 为每个增强子准备独热编码
        one_hot[:, 0] = 1  # 增强子编码为 [1, 0]
    else:
        one_hot = torch.zeros(x.size(0), 2, device=device)  # 为每个增强子准备独热编码
        one_hot[:, 1] = 1  # 增强子编码为 [1, 0]
    
    x = torch.cat([x, one_hot], dim=1)  # 在特征维度上合并类型信息
    
    return x

class EdgeAttentionClassifier(nn.Module):
    """
    边分类器，使用注意力机制加强重要特征，并进行分类。
    """
    def __init__(self, node_features_dim, edge_in_channels, hidden_dim):
        """
        参数:
        node_features_dim (int): 单个节点特征的维度。
        edge_in_channels (int): 边特征的维度。
        hidden_dim (int): 隐藏层的维度。
        """
        super(EdgeAttentionClassifier, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(2 * node_features_dim + edge_in_channels, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Softmax(dim=-1)
        )
        self.classifier = nn.Sequential(
            nn.Linear(2 * node_features_dim + edge_in_channels, hidden_dim),
            nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, combined_features):
        attention_weights = self.attention(combined_features)
        weighted_features = combined_features * attention_weights
        return self.classifier(weighted_features).squeeze()    # 返回 logits

class GATv2Model(nn.Module):
    """
    图注意力网络（GATv2）模型。
    """
    def __init__(self, node_in_channels, edge_in_channels, gat1_out_channels, gat2_out_channels, Edge_hidden_dim, num_heads=1, concat=True, negative_slope=0.2, dropout=0.2, add_self_loops=True, bias=True):
        """
        参数:
        node_in_channels (int): 输入节点特征的维度。
        edge_in_channels (int): 边特征的维度。
        out_channels (int): 每个头部的输出特征维度。
        num_heads (int): 注意力头的数量。
        concat (bool): 如果为True，则在头部间进行连接；否则，取平均。
        negative_slope (float): LeakyReLU中的负斜率。
        dropout (float): 用于注意力系数的dropout率。
        add_self_loops (bool): 是否添加自环。
        bias (bool): 是否添加偏置。
        """
        super(GATv2Model, self).__init__()
        self.gatv2conv1 = GATv2Conv(
            in_channels=node_in_channels,
            out_channels=gat1_out_channels,
            heads=num_heads,
            concat=concat,
            dropout=dropout,
            negative_slope=negative_slope,
            add_self_loops=add_self_loops,
            edge_dim=edge_in_channels,
            bias=bias
        )
        self.gatv2conv2 = GATv2Conv(
            in_channels=gat1_out_channels * num_heads if concat else gat1_out_channels,
            out_channels=gat2_out_channels,
            heads=1,
            concat=False,
            dropout=dropout,
            negative_slope=negative_slope,
            add_self_loops=add_self_loops,
            edge_dim=edge_in_channels,
            bias=bias
        )
        # 边分类器
        self.edge_classifier = EdgeAttentionClassifier(gat2_out_channels, edge_in_channels, hidden_dim=Edge_hidden_dim)

    def forward(self, x, edge_index, edge_attr):
        x = F.dropout(x, p=self.gatv2conv1.dropout, training=self.training)   #通过随机断开一部分神经元的连接来减少过拟合，确保只有在训练模式下应用 dropout。
        # print("Shapes before GATv2Conv1:", x.shape, edge_index.shape, edge_attr.shape)
        x = F.elu(self.gatv2conv1(x, edge_index, edge_attr)) #第一层图神经网络
        # print(x.shape)
        x = F.dropout(x, p=self.gatv2conv2.dropout, training=self.training)   #再次应用 dropout 减少第二层图注意力网络的过拟合风险
        x = self.gatv2conv2(x, edge_index, edge_attr)  #第二层图神经网络
        # print(x.shape)
        # edge_attr = edge_attr.unsqueeze(-1)  # 确保edge_attr维度正确，在 edge_attr 的最后添加一个新维度
        # print(edge_attr.shape, x[edge_index[0]].shape, x[edge_index[1]].shape)
        edge_features = torch.cat([x[edge_index[0]], x[edge_index[1]], edge_attr], dim=-1)   #dim=-1 指在最后一个维度上进行拼接
        # print(edge_features.shape)
        outputs = self.edge_classifier(edge_features)
        if torch.isnan(outputs).any():
            print("NaN detected in model outputs")
        # return torch.sigmoid(self.edge_classifier(edge_features)).squeeze()
        return outputs

    

class DualCNNandGATv2(nn.Module):
    def __init__(self, node_height, node_width, out_channels, kernel_size, cnn_out_channels, gat1_out_channels, gat2_out_channels, Edge_hidden_dim, edge_in_channels, num_heads=1):
        super(DualCNNandGATv2, self).__init__()
        #CNN 在图数据中的操作是针对每一个节点进行的，而不是针对整个图或整个批次,作用于每个节点的特征矩阵
        self.cnn_e = nn.Sequential(
            nn.Conv2d(in_channels=node_height, out_channels=out_channels, kernel_size=(1,kernel_size), padding=(0,(kernel_size - 1) // 2)),  #in_channels是通道,kernel_size=(1,3)高度”维度上的大小为1,在每种特征的21个位置上提取局部模式，而不是跨越不同特征
            nn.LeakyReLU(),
            nn.MaxPool2d(kernel_size=(1,2)),  # 只在宽度维度进行池化
            nn.Flatten(),
            nn.Linear(out_channels * 1 * (node_width // 2) , cnn_out_channels)
        )
        self.cnn_p = nn.Sequential(
            nn.Conv2d(in_channels=node_height, out_channels=out_channels, kernel_size=(1,kernel_size), padding=(0,(kernel_size - 1) // 2)),
            nn.LeakyReLU(),
            nn.MaxPool2d(kernel_size=(1,2)),
            nn.Flatten(),
            nn.Linear(out_channels * 1 * (node_width // 2) , cnn_out_channels)
        )
        self.gatv2 = GATv2Model(node_in_channels=cnn_out_channels+3, edge_in_channels=edge_in_channels, gat1_out_channels=gat1_out_channels, gat2_out_channels=gat2_out_channels, Edge_hidden_dim=Edge_hidden_dim, num_heads=num_heads)

    
    def forward(self, batch_data):
        # 用CNN提取节点特征
        # print(batch_data.x_e.unsqueeze(2).shape)
        x_e = self.cnn_e(batch_data.x_e.unsqueeze(2))  # 从 [B, 7, 21] 变为 [B, 7, 1, 21]
        # print(x_e.shape)
        x_p = self.cnn_p(batch_data.x_p.unsqueeze(2))

        # 拼接节点的度
        degree_e = batch_data.node_degree[:x_e.size(0)]
        degree_p = batch_data.node_degree[x_e.size(0):]
        
        # 拼接节点的类型的编码
        x_e = one_hot_encode(x_e, 'e', x_e.device)
        x_p = one_hot_encode(x_p, 'p', x_p.device)
        
        # 拼接节点的度信息
        x_e = torch.cat([x_e, degree_e], dim=1)
        x_p = torch.cat([x_p, degree_p], dim=1)
        # print(x_e.shape, x_p.shape, batch_data.edge_index.shape, batch_data.edge_attr.shape)

        # 最终节点特征拼接
        x = torch.cat([x_e, x_p], dim=0)  #dim=0上下，dim=1左右
        # 在数据加载或预处理中添加检查
        if torch.isnan(x).any():
            print("NaN detected in input features")

        return self.gatv2(x, batch_data.edge_index, batch_data.edge_attr)

