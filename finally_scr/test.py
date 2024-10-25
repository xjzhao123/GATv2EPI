'''加载最佳超参和模型，测试'''
import torch
from torch_geometric.loader import DataLoader
from GATv2_model import DualCNNandGATv2
from functions import train_model, evaluate_model,load_subgraph_loader, CombinedLoss
from dataset import final_load_datas
import gc  # 导入垃圾收集模块
import os
import json
import pandas as pd
import numpy as np
import random
import time


# # 可以是任意数字，使用相同的数字确保每次运行的初始化相同
# def set_seed(seed):  
#     random.seed(seed)
#     np.random.seed(seed)
#     torch.manual_seed(seed)
#     if torch.cuda.is_available():
#         torch.cuda.manual_seed_all(seed)# 对所有GPU设置种子


# 初始化参数
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
base_dir = '/home/tjzhang03/zxj/deal_data/data_output'
cell_type = ['NHEK', 'IMR90', 'HMEC']  # or ['K562'] for a single cell line
# cell_type = ['K562']
# input_dir = f'{base_dir}/{cell_type}'
# output_dir = f'{input_dir}/random_state_split_dataset_CNN_GAT02'
model_out_dir = f'/home/tjzhang03/zxj/deal_data/result/K562'

geometric_batch_size = 1  #一个batch里面只有一个染色体子图
node_height, node_width = 7, 21  #节点特征的位点，七种特征，21个窗口

# 加载JSON文件
# 定义文件路径
file_path = '/home/tjzhang03/zxj/GATv2EPI/finally_scr/K562_best_hyperparameters.json'
with open(file_path, 'r') as json_file:
    loaded_params = json.load(json_file)

# 解包超参数
out_channels=loaded_params['out_channels']
kernel_size=loaded_params['kernel_size']
cnn_out_channels=loaded_params['cnn_out_channels']
gat1_out_channels=loaded_params['gat1_out_channels']
gat2_out_channels=loaded_params['gat2_out_channels']
edge_hidden_dim=loaded_params['edge_hidden_dim']
num_heads=loaded_params['num_heads']
lr = loaded_params['lr']
weight_decay = loaded_params['weight_decay']
step_size = loaded_params['step_size']
gamma = loaded_params['gamma']



num_epochs = 200
threshold = 0.5

# seed = 0
# set_seed(seed)
aucs_per_state = []
model_save_path = os.path.join(model_out_dir, f'model.pth')
# print(random_state)
#过滤、平衡、划分训练集验证测试、转换数据格式
# train_edges_df, val_edges_df, test_edges_df, train_nodes_dict, val_nodes_dict, test_nodes_dict = create_data(input_dir, cell_type, random_state)
# train_edges_df, val_edges_df, test_edges_df, train_nodes_dict, val_nodes_dict, test_nodes_dict = load_datasets_and_node_features(cell_type, base_dir, random_state)
train_edges_df, test_edges_df, train_nodes_dict, test_nodes_dict = final_load_datas(cell_type, base_dir)

#转换封装成模型需要的输入格式
# train_loader = load_subgraph_loader(train_nodes_dict, train_edges_df, geometric_batch_size)
test_loader = load_subgraph_loader(test_nodes_dict, test_edges_df, geometric_batch_size, is_training=False)


#定义模型
model = DualCNNandGATv2(node_height, node_width, out_channels, kernel_size, cnn_out_channels, gat1_out_channels, 
                        gat2_out_channels, edge_hidden_dim, 2, num_heads).to(device)
model.load_state_dict(torch.load(model_save_path))  #加载模型
criterion = CombinedLoss(weight_bce=0.6)

# 保存每个 epoch 的
train_losses = []
test_losses = []
test_aucs = []

test_loss, test_auc, test_aupr, test_acc = evaluate_model(model, test_loader, criterion, threshold, device)
print(f"Final  test_auc:{test_auc}, test_aupr:{test_aupr}, test_acc:{test_acc}")


# # 循环结束，不再需要模型和数据时
del model  # 删除模型实例
del test_loader # 删除数据变量
gc.collect()  # 显式调用垃圾收集器
if torch.cuda.is_available():
    torch.cuda.empty_cache()  # 清理CUDA缓存


   