import numpy as np
import pandas as pd
import re
import os
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import networkx as nx
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed


def standardize_gene_id(gene_id):
    # 正则表达式去除版本号和其他后缀
    match = re.match(r"(ENSG\d+)", gene_id)
    if match:
        return match.group(1)
    return None

def load_non_gene_ids(file_path):
    # 从文件中加载非基因ID
    with open(file_path, 'r') as file:
        non_gene_ids = set(line.strip() for line in file if line.strip())
    return non_gene_ids


def balance_dataset(input_df, sample_random_state):
    input_df = input_df.copy()
    enhancers_with_positives = input_df[input_df['label'] == 1]['ccRE_ID'].unique()
    filtered_df = input_df[input_df['ccRE_ID'].isin(enhancers_with_positives)]
    balanced_df = pd.DataFrame()
    for ccRE_ID in enhancers_with_positives:
        group = filtered_df[filtered_df['ccRE_ID'] == ccRE_ID]
        pos_samples = group[group['label'] == 1]
        neg_samples = group[group['label'] == 0]
        sampled_neg_samples = neg_samples.sample(n=len(pos_samples), random_state=sample_random_state) if len(neg_samples) >= len(pos_samples) else neg_samples
        balanced_group = pd.concat([pos_samples, sampled_neg_samples])
        balanced_df = pd.concat([balanced_df, balanced_group])
    balanced_df.reset_index(drop=True, inplace=True)
    return balanced_df

def standardize_features(ep_signal, num_exclude):
    columns_to_standardize = ep_signal.columns[:-num_exclude]
    scaler = MinMaxScaler()
    ep_signal[columns_to_standardize] = scaler.fit_transform(ep_signal[columns_to_standardize])
    
    return ep_signal

def convert_to_matrix(df, feature_columns, node_type):
    '''将节点特征向量转换为二维矩阵 (feature,windows)'''
    unique_ids = df['ccRE_ID' if node_type == 'e' else 'Gene_ID'].unique()
    matrices = {}
    feature_bases = set(col.split('_')[0] for col in feature_columns if col.endswith(node_type))

    for node_id in unique_ids:
        node_data = df[df['ccRE_ID' if node_type == 'e' else 'Gene_ID'] == node_id]
        matrix = []
        for feature_base in feature_bases:
            feature_row = [node_data.get(f"{feature_base}_{pos}{node_type}", 0).values[0] if f"{feature_base}_{pos}{node_type}" in node_data.columns else 0 for pos in range(-10, 11)]
            matrix.append(feature_row)
        matrices[node_id] = {
            'matrix': np.array(matrix),
            'chr': node_data['chr'].iloc[0]
        }

    return matrices

def process_ep_signal(ep_signal):
    '''得到节点特征'''
    enhancer_cols = [col for col in ep_signal.columns if col.endswith('e') and col not in ['distance']]
    promoter_cols = [col for col in ep_signal.columns if col.endswith('p') and col not in ['keep']]

    enhancer = ep_signal[enhancer_cols + ['ccRE_ID', 'chr']].drop_duplicates()
    promoter = ep_signal[promoter_cols + ['Gene_ID', 'chr']].drop_duplicates()

    enhancer_features = convert_to_matrix(enhancer, enhancer_cols, 'e')
    promoter_features = convert_to_matrix(promoter, promoter_cols, 'p')
    edges_df = ep_signal[['ccRE_ID', 'Gene_ID', 'label', 'chr', 'distance']].drop_duplicates()

    return enhancer_features, promoter_features, edges_df

def components_to_df(components, edges_df):
    relevant_ids = set()
    for component in components:
        relevant_ids.update(component)
    selected_edges = edges_df[edges_df['ccRE_ID'].isin(relevant_ids) | edges_df['Gene_ID'].isin(relevant_ids)]
    return selected_edges

def process_chromosome(chrom, edges_df, split_random_state):
    chrom_edges_df = edges_df[edges_df['chr'] == chrom]
    G = nx.Graph()
    G.add_nodes_from(chrom_edges_df['ccRE_ID'], bipartite=0)
    G.add_nodes_from(chrom_edges_df['Gene_ID'], bipartite=1)
    G.add_edges_from(zip(chrom_edges_df['ccRE_ID'], chrom_edges_df['Gene_ID']))

    connected_components = list(nx.connected_components(G))
    train_val_components, test_components = train_test_split(connected_components, test_size=0.2, random_state=split_random_state)
    train_components, val_components = train_test_split(train_val_components, test_size=0.25, random_state=split_random_state)  #6:2:2
    
    train_df = components_to_df(train_components, chrom_edges_df)
    val_df = components_to_df(val_components, chrom_edges_df)
    test_df = components_to_df(test_components, chrom_edges_df)
    
    return train_df, val_df, test_df

def split_datasets_by_chromosome(edges_df, output_dir, split_random_state):
    results = Parallel(n_jobs=-1)(delayed(process_chromosome)(chrom, edges_df, split_random_state) for chrom in edges_df['chr'].unique())
    train_dfs, val_dfs, test_dfs = zip(*results)
    
    train_df = pd.concat(train_dfs).drop_duplicates()
    val_df = pd.concat(val_dfs).drop_duplicates()
    test_df = pd.concat(test_dfs).drop_duplicates()
    train_df.to_csv(os.path.join(output_dir, f'{split_random_state}_train_edges.csv'), sep='\t', index=False)
    val_df.to_csv(os.path.join(output_dir, f'{split_random_state}_val_edges.csv'), sep='\t', index=False)
    test_df.to_csv(os.path.join(output_dir, f'{split_random_state}_test_edges.csv'), sep='\t', index=False)

    return train_df, val_df, test_df

def extract_node_ids_from_edges(edges_df):
    ccre_ids = set(edges_df['ccRE_ID'])
    gene_ids = set(edges_df['Gene_ID'])
    node_ids = ccre_ids.union(gene_ids)
    return node_ids


def generate_node_features_for_cnn(train_df, val_df, test_df, enhancer_features, promoter_features, output_dir, random_state):
    datasets = {'train': train_df, 'val': val_df, 'test': test_df}
    result = {}

    for name, dataset in datasets.items():
        enhancer_ids = np.array(list(set(dataset['ccRE_ID'])))
        promoter_ids = np.array(list(set(dataset['Gene_ID'])))

        enhancer_data = {id: enhancer_features[id] for id in enhancer_ids if id in enhancer_features}
        promoter_data = {id: promoter_features[id] for id in promoter_ids if id in promoter_features}

        np.save(os.path.join(output_dir, f'{random_state}_{name}_enhancer_ids.npy'), enhancer_ids)
        np.save(os.path.join(output_dir, f'{random_state}_{name}_enhancer_features.npy'), np.array([data['matrix'] for data in enhancer_data.values()]))
        np.save(os.path.join(output_dir, f'{random_state}_{name}_enhancer_chrs.npy'), np.array([data['chr'] for data in enhancer_data.values()]))
        np.save(os.path.join(output_dir, f'{random_state}_{name}_promoter_ids.npy'), promoter_ids)
        np.save(os.path.join(output_dir, f'{random_state}_{name}_promoter_features.npy'), np.array([data['matrix'] for data in promoter_data.values()]))
        np.save(os.path.join(output_dir, f'{random_state}_{name}_promoter_chrs.npy'), np.array([data['chr'] for data in promoter_data.values()]))

        result[name] = {
            'enhancer': {
                'ids': enhancer_ids,
                'features': np.array([data['matrix'] for data in enhancer_data.values()]),
                'chr': np.array([data['chr'] for data in enhancer_data.values()])
            },
            'promoter': {
                'ids': promoter_ids,
                'features': np.array([data['matrix'] for data in promoter_data.values()]),
                'chr': np.array([data['chr'] for data in promoter_data.values()])
            }
        }

    return result['train'], result['val'], result['test']

# def create_data(input_dir, cell_type, random_state):
#     # 读取数据集
#     output_dir = os.path.join(input_dir, "random_state_split_dataset_CNN_GAT02")
#     os.makedirs(output_dir, exist_ok=True)  
#     # print(output_dir)
#     input_path = os.path.join(input_dir, f"{cell_type}_10windows_epfeature.csv")
#     input_df = pd.read_csv(input_path, sep='\t')
    
#     standardized_df = standardize_features(input_df, num_exclude=10)
#     enhancer_features, promoter_features, edges_df = process_ep_signal(standardized_df)
#     train_edges_df, val_edges_df, test_edges_df = split_datasets_by_chromosome(standardized_df, output_dir, random_state)
#     train_nodes, val_nodes, test_nodes = generate_node_features_for_cnn(train_edges_df, val_edges_df, test_edges_df, enhancer_features, promoter_features, output_dir, random_state)

#     return train_edges_df, val_edges_df, test_edges_df, train_nodes, val_nodes, test_nodes



def load_datasets_and_node_features(output_dir):
    # output_dir = os.path.join(output_dir, "random_state_split_dataset_CNN_GAT02")
    # 加载边数据集
    train_edges_df = pd.read_csv(os.path.join(output_dir, f'train_edges.csv'), sep='\t')
    val_edges_df = pd.read_csv(os.path.join(output_dir, f'val_edges.csv'), sep='\t')
    test_edges_df = pd.read_csv(os.path.join(output_dir, f'test_edges.csv'), sep='\t')

    # 定义辅助函数加载 numpy 保存的数组
    def load_features(dataset_name):
        enhancer_ids = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_ids.npy'), allow_pickle=True)
        enhancer_features = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_features.npy'), allow_pickle=True)
        enhancer_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_chrs.npy'), allow_pickle=True)

        promoter_ids = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_ids.npy'), allow_pickle=True)
        promoter_features = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_features.npy'), allow_pickle=True)
        promoter_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_chrs.npy'), allow_pickle=True)

        return {
            'enhancer': {'ids': enhancer_ids, 'features': enhancer_features, 'chr': enhancer_chrs},
            'promoter': {'ids': promoter_ids, 'features': promoter_features, 'chr': promoter_chrs}
        }

    # 加载节点特征
    train_nodes = load_features('train')
    val_nodes = load_features('val')
    test_nodes = load_features('test')

    return train_edges_df, val_edges_df, test_edges_df, train_nodes, val_nodes, test_nodes


def load_datasets_and_node_features(cell_types, base_dir, save_output_dir):
    combined_train_edges = pd.DataFrame()
    combined_val_edges = pd.DataFrame()
    combined_test_edges = pd.DataFrame()

    # 初始化节点数据存储结构
    combined_train_nodes = {'enhancer': {'ids': [], 'features': [], 'chr': []},
                            'promoter': {'ids': [], 'features': [], 'chr': []}}
    combined_val_nodes = {'enhancer': {'ids': [], 'features': [], 'chr': []},
                          'promoter': {'ids': [], 'features': [], 'chr': []}}
    combined_test_nodes = {'enhancer': {'ids': [], 'features': [], 'chr': []},
                           'promoter': {'ids': [], 'features': [], 'chr': []}}

    for cell_type in cell_types:
        output_dir = os.path.join(base_dir, cell_type, "random_state_split_dataset_CNN_GAT02")
        
        # 加载边数据集
        train_edges_df = pd.read_csv(os.path.join(output_dir, f'train_edges.csv'), sep='\t')
        val_edges_df = pd.read_csv(os.path.join(output_dir, f'val_edges.csv'), sep='\t')
        test_edges_df = pd.read_csv(os.path.join(output_dir, f'test_edges.csv'), sep='\t')
        
        combined_train_edges = pd.concat([combined_train_edges, train_edges_df], ignore_index=True)
        combined_val_edges = pd.concat([combined_val_edges, val_edges_df], ignore_index=True)
        combined_test_edges = pd.concat([combined_test_edges, test_edges_df], ignore_index=True)

        # 加载节点特征
        def append_node_features(features_dict, dataset_name):
            enhancer_ids = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_ids.npy'), allow_pickle=True)
            enhancer_features = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_features.npy'), allow_pickle=True)
            enhancer_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_chrs.npy'), allow_pickle=True)

            promoter_ids = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_ids.npy'), allow_pickle=True)
            promoter_features = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_features.npy'), allow_pickle=True)
            promoter_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_chrs.npy'), allow_pickle=True)

            features_dict['enhancer']['ids'].append(enhancer_ids)
            features_dict['enhancer']['features'].append(enhancer_features)
            features_dict['enhancer']['chr'].append(enhancer_chrs)
            features_dict['promoter']['ids'].append(promoter_ids)
            features_dict['promoter']['features'].append(promoter_features)
            features_dict['promoter']['chr'].append(promoter_chrs)

        append_node_features(combined_train_nodes, 'train')
        append_node_features(combined_val_nodes, 'val')
        append_node_features(combined_test_nodes, 'test')

    # Concatenate all node features along the first axis
    for key in ['enhancer', 'promoter']:
        for sub_key in ['ids', 'features', 'chr']:
            combined_train_nodes[key][sub_key] = np.concatenate(combined_train_nodes[key][sub_key])
            np.save(os.path.join(save_output_dir, f'combined_train_{key}_{sub_key}.npy'), combined_train_nodes[key][sub_key])
            combined_val_nodes[key][sub_key] = np.concatenate(combined_val_nodes[key][sub_key])
            np.save(os.path.join(save_output_dir, f'combined_val_{key}_{sub_key}.npy'), combined_val_nodes[key][sub_key])
            combined_test_nodes[key][sub_key] = np.concatenate(combined_test_nodes[key][sub_key])
            np.save(os.path.join(save_output_dir, f'combined_test_{key}_{sub_key}.npy'), combined_test_nodes[key][sub_key])
    
    # 保存合并后的边数据
    combined_train_edges.to_csv(os.path.join(save_output_dir, f'combined_train_edges.csv'), index=False, sep='\t')
    combined_val_edges.to_csv(os.path.join(save_output_dir, f'combined_val_edges.csv'), index=False, sep='\t')
    combined_test_edges.to_csv(os.path.join(save_output_dir, f'combined_test_edges.csv'), index=False, sep='\t')


    return combined_train_edges, combined_val_edges, combined_test_edges, combined_train_nodes, combined_val_nodes, combined_test_nodes

def final_load_datas(cell_types, base_dir):
    '''合并训练集和验证集'''

    combined_train_edges = pd.DataFrame()
    combined_test_edges = pd.DataFrame()

    # 初始化节点数据存储结构
    combined_train_nodes = {'enhancer': {'ids': [], 'features': [], 'chr': []},
                            'promoter': {'ids': [], 'features': [], 'chr': []}}
    combined_test_nodes = {'enhancer': {'ids': [], 'features': [], 'chr': []},
                           'promoter': {'ids': [], 'features': [], 'chr': []}}

    for cell_type in cell_types:
        output_dir = os.path.join(base_dir, cell_type, "random_state_split_dataset_CNN_GAT02")
        
        # 加载边数据集
        train_edges_df = pd.read_csv(os.path.join(output_dir, f'train_edges.csv'), sep='\t')
        val_edges_df = pd.read_csv(os.path.join(output_dir, f'val_edges.csv'), sep='\t')
        test_edges_df = pd.read_csv(os.path.join(output_dir, f'test_edges.csv'), sep='\t')
        
        # 合并训练集和验证集边数据
        combined_train_edges = pd.concat([combined_train_edges, train_edges_df, val_edges_df], ignore_index=True)
        combined_test_edges = pd.concat([combined_test_edges, test_edges_df], ignore_index=True)

        # 加载节点特征
        def append_node_features(features_dict, dataset_name):
            enhancer_ids = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_ids.npy'), allow_pickle=True)
            enhancer_features = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_features.npy'), allow_pickle=True)
            enhancer_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_enhancer_chrs.npy'), allow_pickle=True)

            promoter_ids = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_ids.npy'), allow_pickle=True)
            promoter_features = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_features.npy'), allow_pickle=True)
            promoter_chrs = np.load(os.path.join(output_dir, f'{dataset_name}_promoter_chrs.npy'), allow_pickle=True)

            features_dict['enhancer']['ids'].append(enhancer_ids)
            features_dict['enhancer']['features'].append(enhancer_features)
            features_dict['enhancer']['chr'].append(enhancer_chrs)
            features_dict['promoter']['ids'].append(promoter_ids)
            features_dict['promoter']['features'].append(promoter_features)
            features_dict['promoter']['chr'].append(promoter_chrs)

        append_node_features(combined_train_nodes, 'train')
        append_node_features(combined_train_nodes, 'val')  
        append_node_features(combined_test_nodes, 'test')

    # Concatenate all node features along the first axis
    for key in ['enhancer', 'promoter']:
        for sub_key in ['ids', 'features', 'chr']:
            combined_train_nodes[key][sub_key] = np.concatenate(combined_train_nodes[key][sub_key])
            combined_test_nodes[key][sub_key] = np.concatenate(combined_test_nodes[key][sub_key])
    
    # # 保存合并后的边数据
    # combined_train_edges.to_csv(os.path.join(save_output_dir, f'{random_state}_combined_train_edges.csv'), index=False, sep='\t')
    # combined_test_edges.to_csv(os.path.join(save_output_dir, f'{random_state}_combined_test_edges.csv'), index=False, sep='\t')

    return combined_train_edges, combined_test_edges, combined_train_nodes, combined_test_nodes

    
