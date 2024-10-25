   
'''划分好的数据集转换为transepi要的格式'''

import pandas as pd
import re
import os
import glob

# # 1. 加载增强子和启动子的注释信息
# def load_annotations(ccRE_file, TSS_file):
#     enhancer_df = pd.read_csv(ccRE_file, sep='\t', header=None,
#                               names=['e_chr', 'enhancer_start', 'enhancer_end', 'drop', 'ccRE_ID', 'ccRE_type'])
#     enhancer_df = enhancer_df.drop(columns=['drop'])

#     tss_df = pd.read_csv(TSS_file, sep='\t', header=None,
#                          names=['p_chr', 'tss', 'tss_dup', 'tss_ID', 'drop', 'strand', 'Gene_ID'])
#     tss_df = tss_df.drop(columns=['tss_dup', 'drop', 'tss_ID'])

#     # print("Loaded TSS columns:", tss_df.columns)  # 打印列名检查
#     print("Unique ccRE_IDs in enhancer data:", enhancer_df['ccRE_ID'].nunique(), "of", len(enhancer_df))
#     print("Unique Gene_IDs in TSS data:", tss_df['Gene_ID'].nunique(), "of", len(tss_df))


#     return enhancer_df, tss_df

# # 2. 映射并处理数据
# def process_data(file, enhancer_df, tss_df):
#     edges_df = pd.read_csv(file, sep='\t')
#     print(edges_df.shape)
#     start_positive_ratio = edges_df['label'].mean()  # 假设正样本标记为print
#     # 映射增强子信息
#     merged_df = edges_df.merge(enhancer_df, how='left', left_on='ccRE_ID', right_on='ccRE_ID')
#     # 映射启动子信息
#     merged_df = merged_df.merge(tss_df, how='left', left_on='Gene_ID', right_on='Gene_ID')

#     # 计算启动子位点
#     merged_df['promoter_start'] = merged_df.apply(lambda row: row['tss'] - 1500 if row['strand'] == '+' else row['tss'] + 500, axis=1)
#     merged_df['promoter_end'] = merged_df.apply(lambda row: row['tss'] + 1500 if row['strand'] == '+' else row['tss'] - 500, axis=1)

#     # 构造输出格式
#     merged_df['e_name'] = merged_df.apply(lambda row: f"{row['e_chr']}:{row['enhancer_start']}-{row['enhancer_end']}|HeLa|{row['ccRE_ID']}", axis=1)
#     merged_df['p_name'] = merged_df.apply(lambda row: f"{row['p_chr']}:{row['promoter_start']}-{row['promoter_end']}|HeLa|{row['Gene_ID']}|{row['tss_ID']}|{row['strand']}", axis=1)

#     final_df = merged_df[['label', 'distance', 'e_chr', 'enhancer_start', 'enhancer_end', 'e_name', 'p_chr', 'promoter_start', 'promoter_end', 'p_name']]

#     # 计算正样本比例
#     over_positive_ratio = final_df['label'].mean()  # 假设正样本标记为1
#     print(f"start_Positive Ratio: {start_positive_ratio:.2f}, over_Positive Ratio: {over_positive_ratio:.2f}")
#     print(final_df.shape)
#     return final_df, over_positive_ratio

# # 3. 主逻辑
# def main(directory, ccRE_file, TSS_file, output_dir, cell_type):
#     enhancer_df, tss_df = load_annotations(ccRE_file, TSS_file)
    
#     for file_path in glob.glob(os.path.join(directory, '*_val_edges.csv')):
#         random_state = os.path.basename(file_path).split('_')[0]
#         processed_data, positive_ratio = process_data(file_path, enhancer_df, tss_df)
#         output_file = os.path.join(output_dir, f'{cell_type}.{random_state}_val.v3.tsv.gz')
#         processed_data.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')
#         # 打印每个文件的正样本比例
#         # print(f"File: {os.path.basename(file_path)}, Positive Ratio: {positive_ratio:.2f}")

# # 路径定义
# cell_type = 'NHEK'
# directory = f'/home/tjzhang03/zxj/deal_data/data_output/{cell_type}/random_state_split_dataset/'
# ccRE_file = '/home/tjzhang03/zxj/BENGI-master/Benchmark/Annotations/hg19-cCREs.bed'
# TSS_file = '/home/tjzhang03/zxj/BENGI-master/Benchmark/Annotations/GENCODEv19-TSSs.bed'
# output_dir = f'/home/tjzhang03/zxj/tansepi/TransEPI-main/data/BENGI_balanced/{cell_type}/'

# # 运行主函数
# main(directory, ccRE_file, TSS_file, output_dir, cell_type)

# 读取原始距离数据并创建字典
def create_distance_dict(file_path):
    df = pd.read_csv(file_path, usecols=['ccRE_ID', 'Gene_ID', 'distance'], sep='\t')
    distance_dict = { (row['ccRE_ID'], row['Gene_ID']): row['distance'] for index, row in df.iterrows()}
    return distance_dict

# 映射并处理数据
def process_TransEPI_data(file, cell_type, distance_dict):
    print("开始process_TransEPI_data")
    edges_df = pd.read_csv(file, sep='\t')
    print(f"Initial data shape: {edges_df.shape}")
    start_positive_ratio = edges_df['label'].mean()  # 计算初始正样本比例

    # 更新距离数据
    edges_df['distance'] = edges_df.apply(lambda row: distance_dict.get((row['ccRE_ID'], row['Gene_ID']), row['distance']), axis=1)
    
    # 构造输出格式
    edges_df['e_name'] = edges_df.apply(lambda row: f"{row['chr']}:{row['enhancer_start']}-{row['enhancer_end']}|{cell_type}|{row['ccRE_ID']}", axis=1)
    edges_df['p_name'] = edges_df.apply(lambda row: f"{row['chr']}:{row['promoter_start']}-{row['promoter_end']}|{cell_type}|{row['Gene_ID']}|{row['tss_ID']}|{row['strand']}", axis=1)

    final_df = edges_df[['label', 'distance', 'chr', 'enhancer_start', 'enhancer_end', 'e_name', 'chr', 'promoter_start', 'promoter_end', 'p_name']]
    final_df.columns = ['label', 'distance', 'e_chr', 'enhancer_start', 'enhancer_end', 'e_name', 'p_chr', 'promoter_start', 'promoter_end', 'p_name']
    
    # 检查第二列distance是否有等于0的值
    if (final_df['distance'] == 0).any():
        print("Warning: There are rows where distance equals 0.")

    # 计算最终正样本比例
    over_positive_ratio = final_df['label'].mean()
    print(f"Start Positive Ratio: {start_positive_ratio:.2f}, Over Positive Ratio: {over_positive_ratio:.2f}")
    print(f"Final data shape: {final_df.shape}")
    return final_df, over_positive_ratio



# def deal_all_file(directory, output_dir, cell_type, distance_dict):
#     '''处理该目录下的所有验证文件'''
#     print("开始")
#     for file_path in glob.glob(os.path.join(directory, '*_test_edges.csv')):
#         random_state = os.path.basename(file_path).split('_')[0]
#         processed_data, positive_ratio = process_TransEPI_data(file_path, cell_type, distance_dict)
#         output_file = os.path.join(output_dir, f'{cell_type}.{random_state}_val.v3.tsv.gz')
#         processed_data.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')


def deal_all_file(directory, output_dir, cell_type, distance_dict):
    print('开始deal_all_file')
    print(directory)
    print(glob.glob(os.path.join(directory, '*_test_edges.csv')))
    # 处理 _test_edges.csv 文件
    for file_path in glob.glob(os.path.join(directory, '*_test_edges.csv')):
        print(f'file_path:{file_path}')
        random_state = os.path.basename(file_path).split('_')[0]
        processed_data, positive_ratio = process_TransEPI_data(file_path, cell_type, distance_dict)
        output_file = os.path.join(output_dir, f'{cell_type}.{random_state}_test.v3.tsv.gz')
        processed_data.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')

    for file_path in glob.glob(os.path.join(directory, '*_train_edges.csv')):
        random_state = os.path.basename(file_path).split('_')[0]
        processed_data, positive_ratio = process_TransEPI_data(file_path, cell_type, distance_dict)
        output_file = os.path.join(output_dir, f'{cell_type}.{random_state}_train.v3.tsv.gz')
        processed_data.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')

    for file_path in glob.glob(os.path.join(directory, '*_val_edges.csv')):
        random_state = os.path.basename(file_path).split('_')[0]
        processed_data, positive_ratio = process_TransEPI_data(file_path, cell_type, distance_dict)
        output_file = os.path.join(output_dir, f'{cell_type}.{random_state}_val.v3.tsv.gz')
        processed_data.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')


# def parse_my_result(file_path):
#     data_list = []
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
    
#     for i in range(0, len(lines), 1):  # 从第二行开始，每隔一行
#         if 'Current random state' in lines[i]:
#             random_state = int(re.search(r'Current random state (\d+)', lines[i]).group(1))
#             max_auc = float(re.search(r'MAX AUC:(\d+\.\d+)', lines[i]).group(1))
#             data_list.append((random_state, max_auc))
    
#     df = pd.DataFrame(data_list, columns=['Random State', 'My Method MAX AUC'])
#     return df

def parse_my_result(file_path):
    data_list = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # 确保lines列表中每两行一组数据处理
    for i in range(0, len(lines), 2):
        # 解析验证集的AUC结果
        val_match = re.search(r'Model saved for random state (\d+) with AUC: (\d+\.\d+)', lines[i])
        if val_match:
            random_state = int(val_match.group(1))
            val_auc = float(val_match.group(2))
        
        # 解析测试集的AUC结果
        test_match = re.search(r'Testing results for random state (\d+): AUC = (\d+\.\d+)', lines[i + 1])
        if test_match:
            # 从逻辑上来说random_state应该和上面解析的一样，这里再次确认
            assert random_state == int(test_match.group(1))
            test_auc = float(test_match.group(2))
        
        # 将结果加入到列表中
        data_list.append((random_state, val_auc, test_auc))
    
    # 将结果列表转换成DataFrame
    df = pd.DataFrame(data_list, columns=['Random State', 'Validation AUC', 'Testing AUC'])
    return df

def parse_random_forest_result(file_path):
    data_list = []
    with open(file_path, 'r') as file:
        for line in file:
            random_state = int(re.search(r'Random State: (\d+)', line).group(1))
            auc_roc = float(re.search(r'AUC-ROC: (\d+\.\d+)', line).group(1))
            data_list.append((random_state, auc_roc))
    
    df = pd.DataFrame(data_list, columns=['Random State', 'Random Forest AUC-ROC'])
    return df

def parse_tranEPI_result(file_path):
    """
    从指定格式的文本文件中解析AUC数据，并计算每个Random State的平均AUC。

    参数:
    file_path (str): 包含数据的文本文件路径。

    返回:
    DataFrame: 包含每个Random State及其平均AUC的DataFrame。
    """
    with open(file_path, 'r') as file:
        data = file.read()

    pattern = r'Running model for random_state (\d+)(.*?)(?=Running model for random_state|$)'
    matches = re.finditer(pattern, data, re.DOTALL)
    
    data_list = []

    for match in matches:
        random_state = int(match.group(1).strip())
        models_data = match.group(2).strip().split('\n')
        auc_values = [float(re.search(r'AUC:\s+(\d+\.\d+)', line).group(1))
                      for line in models_data if 'AUC:' in line]

        if auc_values:  # 检查是否成功解析出AUC值
            for i, auc in enumerate(auc_values):
                data_list.append((random_state, f'Model {i + 1}', auc))
        else:
            print(f"Warning: No AUC values found for random_state {random_state}")

    if not data_list:
        print("No valid data was parsed. Check the file format and the regular expression used.")
        return pd.DataFrame()

    df = pd.DataFrame(data_list, columns=['Random State', 'Model Index', 'AUC'])
    df_pivot = df.pivot(index='Random State', columns='Model Index', values='AUC')
    df_pivot['transEPI_average_AUC'] = df_pivot.mean(axis=1)

    # 重置索引，并仅选择Random State和平均AUC列
    result_df = df_pivot.reset_index()[['Random State', 'transEPI_average_AUC']]

    return result_df

def merge_results(result_dir, analysis_types):
    # 准备文件解析函数的字典
    parsers = {
        'my': parse_my_result,
        'RandomForestClassifier': parse_random_forest_result,
        'transEPI': parse_tranEPI_result
    }

    dataframes = []
    
    # 遍历指定的分析类型，解析对应的文件并收集DataFrame
    for analysis_type in analysis_types:
        file_path = f'{result_dir}/{analysis_type}_result.txt'
        parse_func = parsers.get(analysis_type)
        if parse_func:
            df = parse_func(file_path)
            df.columns = ['Random State'] + [f'{analysis_type}_{col}' for col in df.columns if col != 'Random State']
            dataframes.append(df)

    # 合并所有DataFrame
    final_df = dataframes[0]
    for df in dataframes[1:]:
        final_df = pd.merge(final_df, df, on='Random State')

    # 输出最终的DataFrame
    print(final_df)

    # 保存DataFrame为CSV文件
    final_df.to_csv(f'{result_dir}/final_results.csv', index=False)
    return final_df




def main():

    # # 路径定义
    # directory = '/home/tjzhang03/zxj/deal_data/data_output'
    # cell_type = 'normal'
    # cell_types_normal = ['HMEC', 'IMR90', 'NHEK']  # 定义 normal 对应的细胞系列表
    # input_dir = f'/home/tjzhang03/zxj/deal_data/data_output/{cell_type}/random_state_split_dataset_CNN_GAT02/'
    # output_dir = f'/home/tjzhang03/zxj/tansepi/TransEPI-main/data/BENGI_CNNandGATv2_balanced/{cell_type}/'
    # result_dir = f'/home/tjzhang03/zxj/deal_data/result/{cell_type}'
    # combined_distance_dict = {}  # 初始化一个空字典以存储所有距离数据

    # if cell_type == 'normal':
    #     # 处理每一个 normal 类型的细胞系
    #     for ct in cell_types_normal:
    #         input_file_path = os.path.join(directory, cell_type, f'{ct}_10windows_epfeature.csv')
    #         distance_dict = create_distance_dict(input_file_path)
    #         deal_all_file(input_dir, output_dir, cell_type, distance_dict)

    # else:
    #     # 处理其他类型的细胞系
    #     input_file_path = os.path.join(directory, cell_type, f'{cell_type}_10windows_epfeature.csv')
    #     distance_dict = create_distance_dict(input_file_path)

    # # 生成TransEPI格式的文件
    # deal_all_file(input_dir, output_dir, cell_type, distance_dict)

    # #解析结果文件   #生成文件时，下面都要注释掉
    # analysis_types = ['my', 'transEPI']
    # final_df = merge_results(result_dir, analysis_types)


    # 路径定义
    cell_type = 'NHEK'
    directory = f'/home/tjzhang03/zxj/deal_data/data_output/{cell_type}/random_state_split_dataset_CNN_GAT02/'
    input_file_path = f'/home/tjzhang03/zxj/deal_data/data_output/{cell_type}/{cell_type}_10windows_epfeature.csv'
    output_dir = f'/home/tjzhang03/zxj/tansepi/TransEPI-main/data/BENGI_CNNandGATv2_balanced/{cell_type}/'
    result_dir = f'/home/tjzhang03/zxj/deal_data/result/{cell_type}'

    #得到原始的距离（为经过标准化的）
    distance_dict = create_distance_dict(input_file_path)
    # 生成TransEPI格式的文件
    deal_all_file(directory, output_dir, cell_type, distance_dict)

    # # #解析结果文件   #生成文件时，下面都要注释掉
    # # analysis_types = ['my', 'transEPI']
    # # final_df = merge_results(result_dir, analysis_types)

    


if __name__ == "__main__":
    main()
