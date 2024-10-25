import json
import os
import pickle
import subprocess
from misc_utils import hg19_chromsize, overlap_length
import pybedtools
import pandas as pd
from collections import defaultdict
import numpy as np
from pybedtools import BedTool
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from multiprocessing import Pool, cpu_count
import math
import pandas as pd
from multiprocessing import Pool, cpu_count
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.utils.data as Data
import pandas as pd
import torch.nn.functional as F
import random


# 设置 pybedtools 的临时目录
pybedtools.helpers.set_tempdir('/home/tjzhang03/temp')  

# 检查是否有可用的 GPU
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

# 设置随机数种子
RANDOM_SEED = 42

random.seed(RANDOM_SEED)  # Python随机库
np.random.seed(RANDOM_SEED)  # Numpy随机库
torch.manual_seed(RANDOM_SEED)  # PyTorch随机库

if torch.cuda.is_available():
    torch.cuda.manual_seed(RANDOM_SEED)  # 如果使用GPU，也需要设置



def load_datasets(json_path):
    with open(json_path, 'r') as file:
        data = json.load(file)
    base_location = data['location']
    dataset_paths = [os.path.join(base_location, dataset) for dataset in data['datasets']]
    datasets = {os.path.basename(path): pd.read_csv(path, sep='\t') for path in dataset_paths}
    return datasets

def save_all_dfs(all_dfs, base_dir):
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    file_paths = {}
    for feature_name, df in all_dfs.items():
        file_path = os.path.join(base_dir, f"{feature_name}.csv")
        df.to_csv(file_path)
        file_paths[feature_name] = file_path  # 记录文件路径

    with open(os.path.join(base_dir, "all_dfs_file_paths.json"), 'w') as file:
        json.dump(file_paths, file)

def load_all_dfs(base_dir):
   
    file_paths_path = os.path.join(base_dir, "all_dfs_file_paths.json")
    if not os.path.exists(file_paths_path):
        return None

    with open(file_paths_path, 'r') as file:
        file_paths = json.load(file)

    all_dfs = {}
    for feature_name, file_path in file_paths.items():
        all_dfs[feature_name] = pd.read_csv(file_path, index_col=0)

    return all_dfs

def merge_datasets(datasets):
    combined_dataset = pd.DataFrame()
    for dataset in datasets.values():
        combined_dataset = pd.concat([combined_dataset, dataset], ignore_index=True)
    return combined_dataset


def get_filtered_chrom_sizes():
    """
    获取经过过滤的染色体大小字典。
    移除 'chrY', 'chrM', 'chrMT'。
    """
    chrom_sizes = hg19_chromsize
    keys_to_remove = ['chrY', 'chrM', 'chrMT']
    for key in keys_to_remove:
        chrom_sizes.pop(key, None)  
    return chrom_sizes

def process_bed_files(json_file_path, bin_size):
    with open(json_file_path, 'r') as f:
        data = json.load(f)

    chrom_sizes = get_filtered_chrom_sizes()


    BASE_PATH = data["location"]
    OUTPUT_DIR = data["output_dir"]
    CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

    generated_bed_files = []  

    for cell_type, features in data["celltypes"].items():
        for feature, filename in features.items():
            input_file = os.path.join(BASE_PATH, filename)
            for chrom in CHROMOSOMES:
                chrom_size = chrom_sizes[chrom]
                num_bins = (chrom_size + bin_size - 1) // bin_size
                for i in range(num_bins):
                    start = i * bin_size
                    end = min(start + bin_size, chrom_size)
                    output_file = os.path.join(OUTPUT_DIR, f"{cell_type}_{feature}_{chrom}_{start}_{end}.bed")
                    cmd = f"bigBedToBed -chrom={chrom} -start={start} -end={end} {input_file} {output_file}"
                    subprocess.run(cmd, shell=True)
                    generated_bed_files.append(output_file)

    # 加载.bed文件到数据结构
    signal_dict = load_specific_bed_files(generated_bed_files)

    # 创建BedTool对象
    pybedtools.helpers.set_tempdir('/home/tjzhang03/temp')
    return create_global_bedtool_from_loaded_data(signal_dict)


def load_bed_file_to_data_structure(file_path, data):
    """
    从.bed文件加载数据到嵌套字典结构中
    :param file_path: .bed文件的路径
    :param data: 已存在的数据结构
    """
    # 提取细胞系、修饰特征、染色体和区间信息
    cell_type, feature, chrom, start, end = os.path.basename(file_path).replace('.bed', '').split('_')
    interval = f"{start}_{end}"
    
    # 确保字典中存在相应的键
    if cell_type not in data:
        data[cell_type] = {}
    if feature not in data[cell_type]:
        data[cell_type][feature] = {}
    if chrom not in data[cell_type][feature]:
        data[cell_type][feature][chrom] = {}
    if interval not in data[cell_type][feature][chrom]:
        data[cell_type][feature][chrom][interval] = []

    with open(file_path, 'r') as f:
        for line in f:
            tokens = line.strip().split()
            region_start, region_end, signal = int(tokens[1]), int(tokens[2]), float(tokens[6])
            data[cell_type][feature][chrom][interval].append((region_start, region_end, signal))
    
    return data

def load_specific_bed_files(file_paths):

    data = {}
    for file_path in file_paths:
        try:
            load_bed_file_to_data_structure(file_path, data)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    return data


def create_global_bedtool_from_loaded_data(loaded_data):
  
    global global_signal_bedtool
    bed_entries = [
        (chrom, sig_start, sig_end, cell_type, feature, str(signal_val))  
        for cell_type, cell_data in loaded_data.items()
        for feature, feature_data in cell_data.items()
        for chrom, intervals in feature_data.items()
        for interval, signals in intervals.items()
        for sig_start, sig_end, signal_val in signals
    ]

    # global_signal_bedtool = pybedtools.BedTool(bed_entries)
    return pybedtools.BedTool(bed_entries)  # 返回创建的BedTool对象

def save_bedtool_to_file(bedtool_obj, filename):
    """
    将BedTool对象保存到bed文件。
    """
    bedtool_obj.saveas(filename)

def load_bedtool_from_file(filename):
    """
    从bed文件加载BedTool对象。
    """
    return pybedtools.BedTool(filename)


def calculate_peak_lengths_and_signal_per_feature_and_cell(global_bedtool):
    # 使用defaultdict的双层结构
    cell_feature_data = defaultdict(lambda: defaultdict(list))

    for interval in global_bedtool:
        cell_type = interval.fields[3]
        feature = interval.fields[4]
        peak_length = interval.length
        signal_value = float(interval.fields[5])  # 将信号值转换为浮点数

        cell_feature_data[cell_type][feature].append((peak_length, signal_value))

    return cell_feature_data


def get_peak_and_signal_statistics(cell_feature_data):

    stats = defaultdict(dict)
    for cell_type, feature_data in cell_feature_data.items():
        for feature, data in feature_data.items():
            data_array = np.array(data)
            lengths = data_array[:, 0]
            signal_values = data_array[:, 1]

            stats[cell_type][feature] = {
                'mean_peak_length': int(np.ceil(np.mean(lengths))),  # 使用int()进行转换
                'max_peak_length': int(np.max(lengths)),
                'min_peak_length': int(np.min(lengths)),
                'mean_signal_value': np.mean(signal_values),
                'max_signal_value': np.max(signal_values),
                'min_signal_value': np.min(signal_values)
            }
    return stats


'''全基因组增强子启动子区域，用于后续掩盖'''
def compute_promoter_region(start, end, strand):
    if strand == '+':
        promoter_start = start - 1500
        promoter_end = start + 500
    else:
        promoter_start = end - 500
        promoter_end = end + 1500

    return (min(promoter_start, promoter_end), max(promoter_start, promoter_end))

def get_regions_from_files(enhancer_path, gtf_path):
    unwanted_chroms = ['chrM', 'chrMT', 'chrY']
    enhancers = pybedtools.BedTool(enhancer_path).cut([0, 1, 2])
    enhancers_df = enhancers.to_dataframe()
    enhancers_df = enhancers_df[~enhancers_df['chrom'].isin(unwanted_chroms)]  # 过滤增强子中的不需要的染色体
    enhancers = pybedtools.BedTool.from_dataframe(enhancers_df)

    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, 
                      usecols=[0, 2, 3, 4, 6], names=['chrom', 'tag', 'start', 'end', 'strand'])
    gtf = gtf[~gtf['chrom'].isin(unwanted_chroms)]  #过滤染色体
    gtf = gtf[gtf['tag'] == 'gene']   #过滤['tag'] 为'gene'的

    gtf['start'], gtf['end'] = zip(*gtf.apply(lambda x: compute_promoter_region(x['start'], x['end'], x['strand']), axis=1))

    promoters = pybedtools.BedTool.from_dataframe(gtf[['chrom', 'start', 'end']])    
    regions = enhancers.cat(promoters, postmerge=False).sort()
    return regions

def parse_window_name(window_name):
    index_part = window_name[:-1] 
    type_part = window_name[-1]  
    return (int(index_part), type_part) if index_part else (0, type_part)


def compute_window_positions_test(df, element_type, window_count, window_size, chrom_sizes):
    new_columns = {}

    # 将DataFrame中的相关列转化为NumPy数组，以支持向量化操作。
    start_values = df[f'{element_type}_start'].values
    end_values = df[f'{element_type}_end'].values
    chrom_values = df['chr'].values
    

    # 对每侧的窗口进行迭代。
    for side in ['left', 'right']:
        direction = -1 if side == 'left' else 1

        # 根据窗口数量进行迭代。
        for i in range(1, window_count + 1):
            # 计算新的窗口开始和结束位置。
            if side == 'left':
                new_start = (start_values + direction * i * window_size).clip(min=0)
                new_end = (new_start + window_size).clip(min=0)
            else:
                new_end = (end_values + direction * i * window_size).clip(max=np.array([chrom_sizes[ch] for ch in chrom_values]))
                new_start = (new_end - window_size).clip(max=np.array([chrom_sizes[ch] for ch in chrom_values]))

            # 将新计算的值存储到之前初始化的字典中。
            new_columns[f'{element_type}_{side}_window_{i}_start'] = new_start
            new_columns[f'{element_type}_{side}_window_{i}_end'] = new_end

    # 使用pd.concat一次性将所有新列添加到DataFrame中。
    df = pd.concat([df, pd.DataFrame(new_columns)], axis=1)

    return df


#区域化
def extract_region_info_test(hela_data, window_count):
    all_regions = []

    def build_region_info(row, element_type, side, window_count):
        region_info_list = []
        # 对于增强子和启动子本身，添加自身的区域
        if side == 'self':
            start_col = f'{element_type}_start'
            end_col = f'{element_type}_end'
            side_indicator = 'p' if element_type == 'promoter' else 'e'
            window_name = f"0{side_indicator}"  # 使用0作为编号来表示自身
            region = {
                'chrom': row['chr'],
                'start': row[start_col],
                'end': row[end_col],
                'row_name': row.name,
                'windows_name': window_name
            }
            region_info_list.append(region)
        else:
            for i in range(1, window_count + 1):
                start_col = f'{element_type}_{side}_window_{i}_start'
                end_col = f'{element_type}_{side}_window_{i}_end'
                side_indicator = 'p' if element_type == 'promoter' else 'e'
                window_name = f"{i * (-1 if side == 'left' else 1)}{side_indicator}"

                region = {
                    'chrom': row['chr'],
                    'start': row[start_col],
                    'end': row[end_col],
                    'row_name': row.name,
                    'windows_name': window_name
                }

                region_info_list.append(region)

        return region_info_list

    # 包括自身作为一个区域
    sides = ['left', 'self', 'right']
    for element_type in ['enhancer', 'promoter']:
        for side in sides:
            regions_info = hela_data.apply(build_region_info, element_type=element_type, side=side, window_count=window_count, axis=1)
            all_regions.extend(regions_info)

    all_regions = [item for sublist in all_regions for item in sublist]
    # print("all_regions没问题")
    # 创建 DataFrame
    all_regions_df = pd.DataFrame(all_regions)
    # 排序 DataFrame 行
    all_regions_df = all_regions_df.sort_values(by='windows_name', key=lambda col: col.map(parse_window_name))

    return all_regions_df



# '''可能有点问题'''
# def extract_window_signal_test(signal, all_regions_df):
#     """
#     提取每个窗口的信号值

#     :param signal: 提供的信号数据集
#     :param all_regions_df: 包含所有窗口的DataFrame
#     :return: 一个新的DataFrame，包含窗口的信号值（窗口的各peak的信号值之和）
#     """

#     # 将all_regions_df转换为BedTool对象
#     bed_regions = BedTool.from_dataframe(all_regions_df)

#     # 获取交集
#     intersected = bed_regions.intersect(signal, wa=True, wb=True).to_dataframe()
#     # print(f"Size of intersected DataFrame: {intersected.shape}")

#     # 创建一个字典来跟踪每一行和列名的信号值总和,累积每个窗口的信号值总和。
#     results = {}

#     for row in intersected.itertuples(index=False):
#         idx = int(row[3])  # 获取原本行索引
#         column_name = row[4]  # 获取列名
#         signal_value = float(row[-1])  # 获取信号值

#         # 初始化字典中的条目（如果不存在）
#         if (idx, column_name) not in results:
#             results[(idx, column_name)] = 0  # 字典的键是元组(idx, column_name)

#         # 累加信号值
#         results[(idx, column_name)] += signal_value

#     # 将字典转化为一个DataFrame
#     new_df = pd.DataFrame(list(results.items()), columns=['index_col', 'value'])
#     # 拆分元组列为两列
#     new_df[['idx', 'column_name']] = pd.DataFrame(new_df['index_col'].tolist(), index=new_df.index)
#     # print(f"Size of DataFrame before pivot: {new_df.shape}")

#     # 利用pivot方法进行重构
#     new_df = new_df.pivot(index='idx', columns='column_name', values='value').fillna(0)
#     # print(f"Size of DataFrame before pivot: {new_df.shape}")

#     # 确保所有原始窗口行都包含在最终DataFrame中
#     all_idx = set(all_regions_df['row_name'])  # 使用row_name来保留原始行的索引
#     new_df = new_df.reindex(all_idx).fillna(0)

#     # 按照 unique_windows_names 的顺序重新排序列
#     new_df = new_df.reindex(columns=all_regions_df['windows_name'].unique())

#     return new_df


def extract_window_signal_test(data_df, signal, all_regions_df, mask_regions):
    # 将all_regions_df转换为BedTool对象
    bed_regions = BedTool.from_dataframe(all_regions_df)    #BedTool.from_dataframe()，原本列名被忽略的
    intersected = bed_regions.intersect(signal, wa=True, wb=True).to_dataframe()#前面是bed_regions，后面是bed_regions，列都保留了
    new_columns = [
        'bed_chrom', 'bed_start', 'bed_end', 'row_name', 'windows_name', 
        'sig_chrom', 'sig_start', 'sig_end', 'cell_type', 'feature', 'signal'
    ]
    if intersected.empty:
        # Handle the empty DataFrame, maybe return early or raise an informative error.
        raise ValueError("No intersections found between bed_regions and signal.")
    else:
        intersected.columns = new_columns

    # unique_row_names_in_intersected = intersected['row_name'].nunique()
    # print(f"Unique row names in intersected: {unique_row_names_in_intersected}")
    
    enhancer_promoter_signals = intersected[intersected['windows_name'].isin(["0e", "0p"])]    
    # print(enhancer_promoter_signals.columns)
    flank_regions = intersected[~intersected['windows_name'].isin(["0e", "0p"])]
    flank_bed = BedTool.from_dataframe(flank_regions)      #BedTool.from_dataframe原本的列名被忽略
    filtered_flank = flank_bed.subtract(mask_regions, A=True).to_dataframe()     #A=True -A 模式，有重叠，整个区域删除，不是仅仅剪掉那重叠部分区域
    # print(filtered_flank.columns)
    filtered_flank.columns = new_columns    
    
    # 合并增强子/启动子信号和过滤后的窗口信号
    final_df = pd.concat([enhancer_promoter_signals, filtered_flank], ignore_index=True)   

    # unique_row_names_in_final_df = final_df['row_name'].nunique()
    # print(f"Unique row names in final_df: {unique_row_names_in_final_df}")

    results = {}
    for row in final_df.itertuples(index=False):
        idx = int(row[3])  # 获取原本行索引
        column_name = row[4]  # 获取列名
        signal = float(row[-1])  # 获取信号值

        # 初始化字典中的条目（如果不存在）
        if (idx, column_name) not in results:
            results[(idx, column_name)] = 0

        # 累加信号值
        results[(idx, column_name)] += signal

    full_df = pd.DataFrame(0.0, index=data_df.index, columns=all_regions_df['windows_name'].unique())    #注意得是0.0，因为我的信号值是浮点数
    # 用字典去更新这个full_df
    for (idx, col), value in results.items():
        full_df.at[idx, col] = value


    return full_df

def process_feature(feature_name, cell_type, dataset, window_count, chrom_sizes, global_signal_bedtool, window_size, mask_regions):
    # 这里的 dataset 参数是传入的特定数据集
    data_copy = dataset.copy(deep=True)
    # print('data_copy.shape:',data_copy.shape)

    signal = global_signal_bedtool.filter(lambda x: x[4] == feature_name)

    # 使用先前的函数进行计算，这里只是示例
    data_copy = compute_window_positions_test(data_copy, 'enhancer', window_count, window_size, chrom_sizes)
    data_copy = compute_window_positions_test(data_copy, 'promoter', window_count, window_size, chrom_sizes)
    all_regions_df = extract_region_info_test(data_copy, window_count)
    print("compute_window_positions_test没问题")
    feature_df = extract_window_signal_test(data_copy, signal, all_regions_df, mask_regions)
    print(feature_name,feature_df.shape)

    return feature_name, feature_df

def process_all_features(dataset, cell_type, features, window_count, chrom_sizes, global_signal_bedtool, window_size, mask_regions, output_dir):
    # all_dfs_by_dataset = {}
    # for dataset_name, dataset in datasets.items():
    pool = Pool(processes=cpu_count())

    args = [(feature, cell_type, dataset, window_count, chrom_sizes, global_signal_bedtool, window_size, mask_regions) for feature in features]

    # 使用starmap传递多个参数
    results = pool.starmap(process_feature, args)

    pool.close()
    pool.join()
    
    dfs_by_feature = {feature_name: df for feature_name, df in results}

    updated_dfs = []
    for feature_name, df in dfs_by_feature.items():
        assert not df.isnull().values.any(), f"NaN values detected in feature {feature_name}"
        updated_columns = [f"{feature_name}_{col}" for col in df.columns]
        df.columns = updated_columns
        updated_dfs.append(df)

    # 检查所有DataFrame索引是否一致
    base_index = list(dfs_by_feature.values())[0].index
    for df in dfs_by_feature.values():
        assert all(df.index == base_index), "Inconsistent row index detected!"

    # 将所有的DataFrame按列拼接起来
    final_df = pd.concat(updated_dfs, axis=1)
    final_df["distance"] = dataset["distance"]
    final_df["label"] = dataset["label"]
    final_df["chr"] = dataset["chr"]
    final_df["ccRE_ID"] = dataset["ccRE_ID"]
    final_df["Gene_ID"] = dataset["Gene_ID"]
    final_df["tss_ID"] = dataset["tss_ID"]
    final_df["strand"] = dataset["strand"]
    final_df["enhancer_start"] = dataset["enhancer_start"]
    final_df["enhancer_end"] = dataset["enhancer_end"]
    final_df["promoter_start"] = dataset["promoter_start"]
    final_df["promoter_end"] = dataset["promoter_end"]

    # 保存到文件中
    final_df.to_csv(os.path.join(output_dir, f"{cell_type}_10windows_epfeature.csv"), sep="\t", index=False)

    return final_df





def process_data_for_cell_type_and_(supr_args, cell_type, base_dir):

    print(f'this turn is {cell_type} data')

    # 解包超参数
    BIN_SIZE = supr_args['BIN_SIZE']
    window_count = supr_args['window_count']
    window_size = supr_args['window_size']

    # 路径配置
    genomic_data_dir = f'{base_dir}/genomic_data'
    data_json_dir = f'{base_dir}/data_json'
    process_dir = f'{base_dir}/data_process'
    output_dir = f'{base_dir}/data_output/{cell_type}'
    

    # 其他文件路径配置
    bigbed_json = f'{genomic_data_dir}/json/{cell_type}_bigbed.json'                
    signal_bedtool_file = f'{genomic_data_dir}/genomic_bedtool/{cell_type}_signal_bedtool.bed'       
    enhancer_mask_file = f'{genomic_data_dir}/human_permissive_enhancers_phase_1_and_2.bed'          
    promoter_mask_file = f'{genomic_data_dir}/gencode.v43lift37.annotation.gtf'                     
    MASK_REGIONS_FILE = f'{genomic_data_dir}/genomic_bedtool/mask_EP_regions.bed'                        
    data_json = f'{data_json_dir}/{cell_type}.json'                                    

    datasets = load_datasets(data_json)
    dataset = merge_datasets(datasets)  # 合并所有数据集为一个大的数据集
    
    chrom_sizes = get_filtered_chrom_sizes()

    if os.path.exists(signal_bedtool_file):
        global_signal_bedtool = load_bedtool_from_file(signal_bedtool_file)
    else:
        global_signal_bedtool = process_bed_files(bigbed_json, BIN_SIZE)
        save_bedtool_to_file(global_signal_bedtool, signal_bedtool_file)

    feature_data_per_cell_dict = calculate_peak_lengths_and_signal_per_feature_and_cell(global_signal_bedtool)
    peak_and_signal_statistics = get_peak_and_signal_statistics(feature_data_per_cell_dict)

    if os.path.exists(MASK_REGIONS_FILE):
        mask_regions = load_bedtool_from_file(MASK_REGIONS_FILE)
    else:
        mask_regions = get_regions_from_files(enhancer_mask_file,promoter_mask_file)
        save_bedtool_to_file(mask_regions, MASK_REGIONS_FILE)

    features = set([entry[4] for entry in global_signal_bedtool])

    final_df = process_all_features(dataset, cell_type, features, window_count, chrom_sizes, global_signal_bedtool, window_size, mask_regions, output_dir)

    

def main():
    
    # 定义超参数
    supr_args = {
        'BIN_SIZE': 40000000,  # BIN的大小
        'window_count': 10,
        'window_size': 2000
    }

    base_dir = '/home/tjzhang03/zxj/GATv2EPI'

    # #关于数据选择的参数
    cell_type = "K562"  # 细胞类型,每一个细胞单独处理
    process_data_for_cell_type_and_(supr_args, cell_type, base_dir)


if __name__ == "__main__":
    main()


