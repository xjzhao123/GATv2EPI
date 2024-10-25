
import pandas as pd
import numpy as np
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.data import Batch
from torch.optim import RMSprop
from GATv2_model import DualCNNandGATv2
from functions import train_model, evaluate_model,load_subgraph_loader, CombinedLoss
from dataset import balance_and_split_by_random_state, load_datasets_and_node_features, pre_data_deal, create_data, final_load_datas
import matplotlib.pyplot as plt
import gc  # 导入垃圾收集模块
import random
import os
import json
import time

seed = 42
# 可以是任意数字，使用相同的数字确保每次运行的初始化相同
def set_seed(seed):  
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)# 对所有GPU设置种子


# 初始化参数
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
base_dir = '/home/tjzhang03/zxj/GATv2EPI/data_output'
cell_type = ['NHEK', 'IMR90', 'HMEC']  # or ['K562'] for a single cell line
# input_dir = f'{base_dir}/{cell_type}'
# output_dir = f'{input_dir}/random_state_split_dataset_CNN_GAT02'
model_out_dir = f'/home/tjzhang03/zxj/GATv2EPI/result/normal'
save_output_dir = '/home/tjzhang03/zxj/GATv2EPI/data_output/normal/random_state_split_dataset_CNN_GAT02'

# cell_type = ['K562']
# save_output_dir = '/home/tjzhang03/zxj/deal_data/data_output/K562/random_state_split_dataset_CNN_GAT02'
# model_out_dir = f'/home/tjzhang03/zxj/deal_data/result/K562'

geometric_batch_size = 1  
node_height, node_width = 7, 21  

file_path = '/home/tjzhang03/zxj/GATv2EPI/finally_scr/best_hyperparameters.json'
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

# random_state = 24 
aucs_per_state = []
# train_edges_df, val_edges_df, test_edges_df, train_nodes_dict, val_nodes_dict, test_nodes_dict = load_datasets_and_node_features(cell_type, base_dir, save_output_dir)
model_save_path = os.path.join(model_out_dir, f'model.pth')

train_edges_df, test_edges_df, train_nodes_dict, test_nodes_dict = final_load_datas(cell_type, base_dir)

train_loader = load_subgraph_loader(train_nodes_dict, train_edges_df, geometric_batch_size)
# test_loader = load_subgraph_loader(test_nodes_dict, test_edges_df, geometric_batch_size, is_training=False)

start_time = time.time()  # 记录训练开始时间

#定义模型、优化器等
model = DualCNNandGATv2(node_height, node_width, out_channels, kernel_size, cnn_out_channels, gat1_out_channels, 
                        gat2_out_channels, edge_hidden_dim, 2, num_heads).to(device)
optimizer = RMSprop(model.parameters(), lr=lr, weight_decay=weight_decay)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=gamma)
criterion = CombinedLoss(weight_bce=0.6)

# 保存每个 epoch 的
best_test_auc = float('-inf')
last_train_loss = float('inf')
# last_val_loss = float('inf')
best_model = None

train_losses = []
test_losses = []
test_aucs = []

for epoch in range(num_epochs):
    # 训练和验证模型的代码
    train_loss, train_acc = train_model(model, train_loader, criterion, optimizer, threshold, device)
    train_losses.append(train_loss)
    scheduler.step()  # 更新学习率

# 训练循环结束后
total_training_time = time.time() - start_time
print("Total training time: {:.2f} seconds".format(total_training_time))


# 保存模型
torch.save(model.state_dict(), model_save_path)

# # 循环结束，不再需要模型和数据时
del model  # 删除模型实例
del train_loader# 删除数据变量
gc.collect()  # 显式调用垃圾收集器
if torch.cuda.is_available():
    torch.cuda.empty_cache()  # 清理CUDA缓存




# '''K562'''
# import pandas as pd
# import numpy as np
# import torch
# from torch_geometric.loader import DataLoader
# from torch_geometric.data import Batch
# from torch.optim import RMSprop
# from GATv2_model import DualCNNandGATv2
# from functions import train_model, evaluate_model,load_subgraph_loader, CombinedLoss
# from dataset import load_datasets_and_node_features
# import matplotlib.pyplot as plt
# import gc  # 导入垃圾收集模块
# import random
# import os
# import time


# seed = 42  # 可以是任意数字，使用相同的数字确保每次运行的初始化相同
# torch.manual_seed(seed)
# np.random.seed(seed)
# random.seed(seed)
# if torch.cuda.is_available():
#     torch.cuda.manual_seed_all(seed)  # 对所有GPU设置种子


# # 初始化参数
# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
# base_dir = '/home/tjzhang03/zxj/GATv2EPI/data_output'
# cell_type = ['K562']
# save_output_dir = '/home/tjzhang03/zxj/GATv2EPI/data_output/K562/random_state_split_dataset_CNN_GAT02'
# model_out_dir = f'/home/tjzhang03/zxj/GATv2EPI/result/K562'

# geometric_batch_size = 1  #一个batch里面只有一个染色体子图
# node_height, node_width = 7, 21  #节点特征的位点，七种特征，21个窗口

# out_channels = 32 
# kernel_size = 3
# cnn_out_channels = 64  #节点通过CNN后的特征维度 64
# gat1_out_channels = 32  #第一层GAT输出维度
# gat2_out_channels = 16  #第二层GAT输出节点维度
# Edge_hidden_dim = 16
# num_heads=6
# lr=0.0001
# weight_decay=1e-4
# step_size=30
# gamma=0.8
# num_epochs = 200
# threshold = 0.5


# random_states = [21] # 示例：遍历5个不同的random_state

# aucs_per_state = []
# model_save_path = os.path.join(model_out_dir, f'model.pth')
#     # print(random_state)
#     #过滤、平衡、划分训练集验证测试、转换数据格式
# train_edges_df, val_edges_df, test_edges_df, train_nodes_dict, val_nodes_dict, test_nodes_dict = load_datasets_and_node_features(cell_type, base_dir, save_output_dir)

# #转换封装成模型需要的输入格式
# train_loader = load_subgraph_loader(train_nodes_dict, train_edges_df, geometric_batch_size)
# val_loader = load_subgraph_loader(val_nodes_dict, val_edges_df, geometric_batch_size, is_training=False) #val_loader是字典
# test_loader = load_subgraph_loader(test_nodes_dict, test_edges_df, geometric_batch_size, is_training=False)

# start_time = time.time()  # 记录训练开始时间
# #定义模型、优化器等
# model = DualCNNandGATv2(node_height, node_width, out_channels, kernel_size, cnn_out_channels, gat1_out_channels, 
#                         gat2_out_channels, Edge_hidden_dim, 2, num_heads).to(device)
# optimizer = RMSprop(model.parameters(), lr=lr, weight_decay=weight_decay)
# scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=gamma)
# criterion = CombinedLoss(weight_bce=0.6)

# # 保存每个 epoch 的
# best_val_auc = float('-inf')
# last_train_loss = float('inf')
# last_val_loss = float('inf')
# best_model = None

# train_losses = []
# val_losses = []
# val_aucs = []

# for epoch in range(num_epochs):
#     # 训练和验证模型的代码
#     train_loss, train_acc = train_model(model, train_loader, criterion, optimizer, threshold, device)
#     train_losses.append(train_loss)
#     val_loss, val_auc, _, _ = evaluate_model(model, val_loader, criterion, threshold, device)
#     val_losses.append(val_loss)
#     val_aucs.append(val_auc)
#     scheduler.step()  # 更新学习率

#     if val_auc > best_val_auc:
#         best_val_auc = val_auc
#         best_model = model.state_dict()  # Save the model state

# # 训练循环结束后
# total_training_time = time.time() - start_time
# print("Total training time: {:.2f} seconds".format(total_training_time))
# #保存AUC最高那个模型模型
# # model_save_path = os.path.join(model_out_dir, f'model_{random_state}.pth')
# torch.save(best_model, model_save_path)
# print(f"Model saved with AUC: {best_val_auc}")

# # Load best model for testing 测试
# model.load_state_dict(best_model)
# test_loss, test_auc, test_aupr, test_acc = evaluate_model(model, test_loader, criterion, threshold, device)
# print(f"Final test_auc:{test_auc}, test_aupr:{test_aupr}, test_acc:{test_acc}")


# # # 循环结束，不再需要模型和数据时
# del model  # 删除模型实例
# del train_loader, val_loader, test_loader # 删除数据变量
# gc.collect()  # 显式调用垃圾收集器
# if torch.cuda.is_available():
#     torch.cuda.empty_cache()  # 清理CUDA缓存
