


'''EG2EP'''
# import pandas as pd
# import os

# # 指定数据集文件所在的目录和结果保存目录
# # data_directory = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/Remove-Ambiguous-Pairs.Natural-Ratio/conbine_data/'
# data_directory = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/All-Pairs.Natural-Ratio/'
# # data_directory = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/Remove-Ambiguous-Pairs.Natural-Ratio/'
# output_directory = '/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio'

# # 读取增强子和TSS注释文件，并创建对应字典
# df_cCREs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/hg19-cCREs.bed', sep='\t', usecols=[0, 1, 2, 4, 5], header=None, names=['chr', 'enhancer_start', 'enhancer_end', 'ccRE_ID', 'ccRE_type'])
# cCREs_dict = dict(zip(df_cCREs["ccRE_ID"], df_cCREs[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]].values))

# df_TSSs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/GENCODEv19-TSSs.bed', sep='\t', usecols=[0, 1, 3, 5, 6], header=None, names=['chr', 'tss', 'tss_ID', 'strand', 'Gene_ID'])
# p_dict = dict(zip(df_TSSs["Gene_ID"], df_TSSs[["chr", "tss", "tss_ID", "strand"]].values))

# # 函数计算启动子区域(区分正负链，TSS上游1500bp到下游500bp)
# def compute_promoter_region(tss, strand):
#     if strand == '+':
#         return tss - 1500, tss + 500
#     else:
#         return tss - 500, tss + 1500

# # 遍历数据文件
# for data_file in os.listdir(data_directory):
#     if data_file.endswith('.txt'):
#         file_path = os.path.join(data_directory, data_file)      #文件路径和文件名
#         df1 = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3], header=None, names=['ccRE_ID', 'Gene_ID', 'label', 'cv'])       #读取EG数据集

#         # 向量化操作（在整个数组或数据结构上一次性执行操作）获取增强子和启动子数据
#         df1["enhancer_data"] = df1["ccRE_ID"].map(cCREs_dict)     #map()用于将一个映射（例如字典）应用到 Series 中的每个元素
#         df1["promoter_data"] = df1["Gene_ID"].map(p_dict)
#         df1.dropna(subset=["enhancer_data", "promoter_data"], inplace=True)      #删除含有缺失值（NaN）的行

#         # 拆分增强子和启动子数据
#         df1[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]] = pd.DataFrame(df1["enhancer_data"].tolist(), index=df1.index)
#         df1[["tss_chr", "tss", "tss_ID", "strand"]] = pd.DataFrame(df1["promoter_data"].tolist(), index=df1.index)
#         df1.drop(columns=["enhancer_data", "promoter_data"], inplace=True)

#         # 计算启动子区域
#         df1["promoter_start"], df1["promoter_end"] = zip(*df1.apply(lambda x: compute_promoter_region(x["tss"], x["strand"]), axis=1))
        
#         # 计算距离
#         # df1["distance"] = df1.apply(lambda row: row['promoter_start'] - row['enhancer_end'] if row['enhancer_start'] < row['promoter_start'] else row['enhancer_start'] - row['promoter_end'], axis=1)
#         df1["distance"] = df1.apply(lambda row: abs((row['enhancer_start'] + row['enhancer_end'])/2 - row['tss']), axis=1)
#         # 选择需要的数据并保存
#         df_out = df1[df1["ccRE_type"] == "Enhancer-like"][["chr", "enhancer_start", "enhancer_end", "promoter_start", "promoter_end", "distance", "label", "cv", "strand", "ccRE_ID", "Gene_ID", "tss_ID", "ccRE_type"]]
#         # df_out.rename(columns={'tag': 'label'}, inplace=True)  # 重命名列名 tag 为 label
#         ep_output_file = os.path.join(output_directory, f'{os.path.splitext(data_file)[0]}_ep_distanceisWindows.txt')
#         df_out.to_csv(ep_output_file, sep="\t", index=None, header=True)#包含列名，以后就可以直接读取，不用再给列命名了



# import pandas as pd
# import os



# # 选择基于链方向的TSS
# def select_best_tss(tss_list):
#     if not tss_list:
#         return None
#     # 按照TSS位置排序
#     tss_list.sort(key=lambda x: x[1])
#     if tss_list[0][3] == '+':  # 正链，选择第一个TSS
#         return tss_list[0]
#     else:  # 负链，选择最后一个TSS
#         return tss_list[-1]


# # 函数计算启动子区域(区分正负链，TSS上游1500bp到下游500bp)
# def compute_promoter_region(tss, strand):
#     if strand == '+':
#         return tss - 1500, tss + 500
#     else:
#         return tss - 500, tss + 1500

# data_directory = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/All-Pairs.Natural-Ratio/'
# output_directory = '/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio'

# # 读取数据
# df_cCREs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/hg19-cCREs.bed', sep='\t', usecols=[0, 1, 2, 4, 5], header=None, names=['chr', 'enhancer_start', 'enhancer_end', 'ccRE_ID', 'ccRE_type'])
# cCREs_dict = dict(zip(df_cCREs["ccRE_ID"], df_cCREs[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]].values))

# df_TSSs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/GENCODEv19-TSSs.bed', sep='\t', usecols=[0, 1, 3, 5, 6], header=None, names=['chr', 'tss', 'tss_ID', 'strand', 'Gene_ID'])
# p_dict = dict(zip(df_TSSs["Gene_ID"], df_TSSs[["chr", "tss", "tss_ID", "strand"]].values))


# total_samples = 0
# total_enhancer_length = 0

# # 遍历数据文件
# for data_file in os.listdir(data_directory):
#     if data_file.endswith('.txt'):
#         file_path = os.path.join(data_directory, data_file)
#         df1 = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3], header=None, names=['ccRE_ID', 'Gene_ID', 'label', 'cv'])

#         df1["enhancer_data"] = df1["ccRE_ID"].map(cCREs_dict)
#         df1["promoter_data"] = df1["Gene_ID"].map(p_dict)
#         df1.dropna(subset=["enhancer_data", "promoter_data"], inplace=True)

#         df1[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]] = pd.DataFrame(df1["enhancer_data"].tolist(), index=df1.index)
#         df1[["tss_chr", "tss", "tss_ID", "strand"]] = pd.DataFrame(df1["promoter_data"].tolist(), index=df1.index)
#         df1.drop(columns=["enhancer_data", "promoter_data"], inplace=True)

#         df1["enhancer_length"] = df1["enhancer_end"] - df1["enhancer_start"]
#         df1["promoter_start"], df1["promoter_end"] = zip(*df1.apply(lambda x: compute_promoter_region(x["tss"], x["strand"]), axis=1))

#         # 计算距离
#         df1["distance"] = df1.apply(lambda row: abs((row['enhancer_start'] + row['enhancer_end'])/2 - row['tss']), axis=1)

#         df_out = df1[df1["ccRE_type"] == "Enhancer-like"][["chr", "enhancer_start", "enhancer_end", "promoter_start", "promoter_end", "distance", "label", "cv", "strand", "ccRE_ID", "Gene_ID", "tss_ID", "ccRE_type", "enhancer_length"]]

#         # 计算当前文件的平均增强子长度
#         current_avg_enhancer_length = df_out['enhancer_length'].mean()
#         print(f'{data_file}: Sample size={len(df_out)}, Average Enhancer Length={current_avg_enhancer_length}')
        
#         # 累加总样本量和总长度
#         total_samples += len(df_out)
#         total_enhancer_length += df_out['enhancer_length'].sum()

#         # 删除不再需要的列
#         df_out.drop(columns=["enhancer_length"], inplace=True)
#         ep_output_file = os.path.join(output_directory, f'{os.path.splitext(data_file)[0]}_ep_distanceisWindows.txt')
#         df_out.to_csv(ep_output_file, sep="\t", index=None, header=True) #包含列名，以后就可以直接读取，不用再给列命名了

# # 计算整体平均增强子长度和启动子长度
# overall_avg_enhancer_length = total_enhancer_length / total_samples

# print(f'Overall Average Enhancer Length: {overall_avg_enhancer_length}')

# # 根据你的说明，增强子和启动子的平均长度
# average_length = (overall_avg_enhancer_length + 1500) / 2
# print(f'Combined Average Length of Enhancers and Promoters: {average_length}')





# import pandas as pd
# import os

# # 选择基于链方向的TSS  一个基因好多个TSS
# def select_best_tss(tss_list):
#     if not tss_list:
#         return None
#     # 按照TSS位置排序
#     tss_list.sort(key=lambda x: x[1])
#     if tss_list[0][3] == '+':  # 正链，选择第一个TSS
#         return tss_list[0]
#     else:  # 负链，选择最后一个TSS
#         return tss_list[-1]


# # 计算启动子区域（考虑正负链）
# def compute_promoter_region(tss, strand):
#     if strand == '+':
#         return tss - 1500, tss + 500
#     else:
#         return tss - 500, tss + 1500

# # 读取数据
# df_cCREs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/hg19-cCREs.bed', sep='\t', usecols=[0, 1, 2, 4, 5], header=None, names=['chr', 'enhancer_start', 'enhancer_end', 'ccRE_ID', 'ccRE_type'])
# df_TSSs = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_input/BENGI/GENCODEv19-TSSs.bed', sep='\t', usecols=[0, 1, 3, 5, 6], header=None, names=['chr', 'tss', 'tss_ID', 'strand', 'Gene_ID'])

# # 创建从Gene_ID到所有TSS的映射
# p_dict = {}
# for _, row in df_TSSs.iterrows():
#     gene_id = row['Gene_ID']
#     tss_info = (row['chr'], row['tss'], row['tss_ID'], row['strand'])
#     if gene_id not in p_dict:
#         p_dict[gene_id] = []
#     p_dict[gene_id].append(tss_info)
# #映射
# cCREs_dict = dict(zip(df_cCREs["ccRE_ID"], df_cCREs[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]].values))
# for gene_id in p_dict:
#     p_dict[gene_id] = select_best_tss(p_dict[gene_id])

# data_directory = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/All-Pairs.Natural-Ratio/'
# output_directory = '/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio'

# total_samples = 0
# total_enhancer_length = 0


# # 遍历数据文件
# for data_file in os.listdir(data_directory):
#     if data_file.endswith('.txt'):
#         file_path = os.path.join(data_directory, data_file)
#         df1 = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3], header=None, names=['ccRE_ID', 'Gene_ID', 'label', 'cv'])

#         df1["enhancer_data"] = df1["ccRE_ID"].map(cCREs_dict)
#         df1["promoter_data"] = df1["Gene_ID"].map(p_dict)
#         df1.dropna(subset=["enhancer_data", "promoter_data"], inplace=True)

#         df1[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]] = pd.DataFrame(df1["enhancer_data"].tolist(), index=df1.index)
#         df1[["tss_chr", "tss", "tss_ID", "strand"]] = pd.DataFrame(df1["promoter_data"].tolist(), index=df1.index)
#         df1.drop(columns=["enhancer_data", "promoter_data"], inplace=True)

#         df1["enhancer_length"] = df1["enhancer_end"] - df1["enhancer_start"]
#         df1["promoter_start"], df1["promoter_end"] = zip(*df1.apply(lambda x: compute_promoter_region(x["tss"], x["strand"]), axis=1))

#         # 计算距离
#         df1["distance"] = df1.apply(lambda row: abs((row['enhancer_start'] + row['enhancer_end'])/2 - row['tss']), axis=1)

#         df_out = df1[df1["ccRE_type"] == "Enhancer-like"][["chr", "enhancer_start", "enhancer_end", "promoter_start", "promoter_end", "distance", "label", "cv", "strand", "ccRE_ID", "Gene_ID", "tss_ID", "ccRE_type", "enhancer_length"]]

#         # 计算当前文件的平均增强子长度
#         current_avg_enhancer_length = df_out['enhancer_length'].mean()
#         print(f'{data_file}: Sample size={len(df_out)}, Average Enhancer Length={current_avg_enhancer_length}')
        
#         # 累加总样本量和总长度
#         total_samples += len(df_out)
#         total_enhancer_length += df_out['enhancer_length'].sum()

#         # 删除不再需要的列
#         df_out.drop(columns=["enhancer_length"], inplace=True)
#         ep_output_file = os.path.join(output_directory, f'{os.path.splitext(data_file)[0]}_ep_distanceisWindows.txt')
#         df_out.to_csv(ep_output_file, sep="\t", index=None, header=True) #包含列名，以后就可以直接读取，不用再给列命名了


# # 计算整体平均增强子长度
# overall_avg_enhancer_length = total_enhancer_length / total_samples if total_samples > 0 else 0
# print(f'Overall Average Enhancer Length: {overall_avg_enhancer_length}')



'''EP2EG'''
import pandas as pd
import numpy as np
import os
import re
from misc_utils import hg19_chromsize

# 选择距离增强子中点最近的TSS,一个基因好多个TSS
def select_best_tss(tss_list, enhancer_midpoint):
    if not tss_list:
        return None
    # 选择距离增强子中点最近的TSS
    return min(tss_list, key=lambda x: abs(x[1] - enhancer_midpoint))

def compute_promoter_region(tss, strand, chr_size):
    if strand == '+':
        start = max(1, tss - 1500)
        end = min(tss + 500, chr_size)
    else:
        start = max(1, tss - 500)
        end = min(tss + 1500, chr_size)
    return start, end

def standardize_id(identifier):
    # 正则表达式去除版本号和其他后缀
    match = re.match(r"(ENS[G|T]\d+)", identifier)
    if match:
        return match.group(1)
    return None

def parse_attributes(attribute_str):
    # 解析属性字段
    attributes = {}
    for attr in attribute_str.split(';'):
        if attr.strip():
            key_value = attr.strip().split(' ')
            if len(key_value) == 2:
                key, value = key_value
                attributes[key.strip()] = value.strip('"')
    return attributes

def load_and_process_gtf(gtf_path, chr_size_dict):
    '''从全基因组文件初步筛选TSS'''
    # 加载 GTF 文件
    gtf_data = pd.read_csv(gtf_path, comment='#', sep='\t', header=None, 
                           names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # 筛选出转录本
    transcripts = gtf_data[gtf_data['feature'] == 'transcript']
    transcripts['attributes'] = transcripts['attribute'].apply(parse_attributes)
    transcripts = transcripts[(transcripts['attributes'].apply(lambda x: x.get('transcript_type') in ['protein_coding', 'lincRNA', 'antisense_RNA'])) &
                              (transcripts['attributes'].apply(lambda x: int(x.get('level', 3)) <= 2))]

    # 提取 TSS 并计算启动子区域
    transcripts['TSS'] = transcripts.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)
    transcripts['chr_size'] = transcripts['chr'].map(chr_size_dict)
    transcripts[['promoter_start', 'promoter_end']] = transcripts.apply(lambda x: compute_promoter_region(x['TSS'], x['strand'], x['chr_size']), axis=1, result_type='expand')

    # 标准化 gene_id 和 transcript_id
    transcripts['Gene_ID'] = transcripts['attributes'].apply(lambda x: standardize_id(x.get('gene_id')))
    transcripts['tss_ID'] = transcripts['attributes'].apply(lambda x: standardize_id(x.get('transcript_id')))

    # 选择需要的列
    final_transcripts = transcripts[['chr', 'promoter_start', 'promoter_end', 'strand', 'TSS', 'Gene_ID', 'tss_ID']]
    return final_transcripts



def generate_negative_samples(enhancer_row, transcripts_df, positive_genes, sample_counts):
    '''每个增强子，筛选TSS'''
    # 筛选条件
    mask = (
        (transcripts_df['chr'] == enhancer_row['chr']) &
        (~transcripts_df['Gene_ID'].isin(positive_genes)) &
        (42000 <= np.abs(transcripts_df['TSS'] - enhancer_row['enhancer_midpoint'])) &
        (np.abs(transcripts_df['TSS'] - enhancer_row['enhancer_midpoint']) <= 500000)
    )
    valid_transcripts = transcripts_df[mask]
    num_pos_samples = sample_counts.get(enhancer_row['ccRE_ID'], 0)

    # 检查可用的负样本数量并调整请求的样本数量
    num_samples_to_return = min(len(valid_transcripts), num_pos_samples)
    # if num_samples_to_return < num_pos_samples:
        # print(f"Not enough negative samples for CCRE_ID {enhancer_row['ccRE_ID']}. Expected {num_pos_samples}, got {num_samples_to_return}.")
    
    # 如果有足够的负样本候选，随机选择与正样本数量相同的负样本
    if len(valid_transcripts) > 0:
        samples = valid_transcripts.sample(n=num_samples_to_return)
        samples['enhancer_start'] = enhancer_row['enhancer_start']
        samples['enhancer_end'] = enhancer_row['enhancer_end']
        samples['label'] = 0  # 设置负样本标签为0
        samples['distance'] = np.abs(samples['TSS'] - enhancer_row['enhancer_midpoint'])
        samples['ccRE_ID'] = enhancer_row['ccRE_ID']
        return samples[['chr', 'enhancer_start', 'enhancer_end', 'promoter_start', 'promoter_end', 'distance', 'label', 'strand', 'ccRE_ID', 'Gene_ID', 'tss_ID']]
    else:
        return pd.DataFrame()
    
def generate_negative_samples_for_all(positive_enhancers, transcripts_df, positive_genes, sample_counts):
    '''为所有增强子生成负样本'''
    negative_samples = pd.DataFrame()
    generated_samples = set()  ##避免对同一 CCRE_ID 多次生成负样本
    # 遍历每个正样本增强子
    for index, enhancer in positive_enhancers.iterrows():
        if enhancer['ccRE_ID'] not in generated_samples:
            #调用生成负样本的函数
            neg_samples = generate_negative_samples(enhancer, transcripts_df, positive_genes, sample_counts)
            negative_samples = pd.concat([negative_samples, neg_samples], ignore_index=True)
            generated_samples.add(enhancer['ccRE_ID'])        
    return negative_samples

def process_and_save_samples(positive_examples, transcripts_df, ep_output_file):
    '''得到负样本，平衡的完整数据集'''
    # 加载数据
    # positive_examples = pd.read_csv(input_file_path, sep='\t')
    positive_examples['Gene_ID'] = positive_examples['Gene_ID'].apply(standardize_id)
    positive_examples['tss_ID'] = positive_examples['tss_ID'].apply(standardize_id)

    # 筛选出距离范围内的样本
    pos_samples_df = positive_examples[(positive_examples['distance'] >= 42000) & (positive_examples['distance'] <= 500000)]
    pos_samples_df = pos_samples_df[['chr', 'enhancer_start', 'enhancer_end', 'promoter_start', 'promoter_end', 'distance', 'label', 'strand', 'ccRE_ID', 'Gene_ID', 'tss_ID', 'enhancer_midpoint']]

    # 预计算每个增强子的正样本数量
    sample_counts = pos_samples_df['ccRE_ID'].value_counts()

    # 获取增强子信息，避免后续重复的drop_duplicates操作
    positive_enhancers = pos_samples_df.drop_duplicates()
    # 使用set存储Gene_ID，用于快速查找和检查
    positive_genes = set(positive_enhancers['Gene_ID'])  

    # 生成负样本
    negative_samples = generate_negative_samples_for_all(positive_enhancers, transcripts_df, positive_genes, sample_counts)
    print(f"Number of Positive Samples: {pos_samples_df.shape[0]}")
    print(f"Number of Negative Samples: {negative_samples.shape[0]}")
    # 合并正负样本并按照ccRE_ID进行分组，确保每组内正样本紧跟其负样本
    all_samples = pd.concat([pos_samples_df, negative_samples], ignore_index=True)
    all_samples_sorted = all_samples.sort_values(by=['ccRE_ID', 'label'], ascending=[True, False])

    # 保存DataFrame到CSV
    all_samples_sorted.to_csv(ep_output_file, sep="\t", index=None, header=True)


def main():

    #定义路径
    data_directory = '/home/tjzhang03/zxj/GATv2EPI/data_input/BENGI/All-Pairs.Natural-Ratio/'
    # output_directory = '/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio'
    output_dir = '/home/tjzhang03/zxj/GATv2EPI/data_process/create_data'

    #从全基因组注释文件annotation.gtf筛选TSS，为生成负样本做准备
    gtf_path = '/home/tjzhang03/zxj/GATv2EPI/genomic_data/gencode.v43lift37.annotation.gtf'
    transcripts_df = load_and_process_gtf(gtf_path, hg19_chromsize)

    # 读取数据集注释数据
    df_cCREs = pd.read_csv('/home/tjzhang03/zxj/GATv2EPI/data_input/BENGI/hg19-cCREs.bed', sep='\t', usecols=[0, 1, 2, 4, 5], header=None, names=['chr', 'enhancer_start', 'enhancer_end', 'ccRE_ID', 'ccRE_type'])
    df_TSSs = pd.read_csv('/home/tjzhang03/zxj/GATv2EPI/data_input/BENGI/GENCODEv19-TSSs.bed', sep='\t', usecols=[0, 1, 3, 5, 6], header=None, names=['chr', 'tss', 'tss_ID', 'strand', 'Gene_ID'])

    # 创建从Gene_ID到所有TSS的映射
    p_dict = {}
    for _, row in df_TSSs.iterrows():
        gene_id = row['Gene_ID']
        tss_info = (row['chr'], row['tss'], row['tss_ID'], row['strand'])
        if gene_id not in p_dict:
            p_dict[gene_id] = []
        p_dict[gene_id].append(tss_info)
    #映射
    cCREs_dict = dict(zip(df_cCREs["ccRE_ID"], df_cCREs[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]].values))

    total_samples = 0
    total_enhancer_length = 0

    # 遍历数据文件
    for data_file in os.listdir(data_directory):
        if data_file.endswith('.txt'):
            file_path = os.path.join(data_directory, data_file)
            df1 = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3], header=None, names=['ccRE_ID', 'Gene_ID', 'label', 'cv'])
            #获取增强子信息
            df1["enhancer_data"] = df1["ccRE_ID"].map(cCREs_dict)
            df1.dropna(subset=["enhancer_data"], inplace=True)
            df1[["chr", "enhancer_start", "enhancer_end", "ccRE_type"]] = pd.DataFrame(df1["enhancer_data"].tolist(), index=df1.index)
            # 计算增强子中点
            df1["enhancer_midpoint"] = (df1["enhancer_start"] + df1["enhancer_end"]) / 2
            # 根据增强子中点选择最近的TSS
            df1["promoter_data"] = df1.apply(lambda x: select_best_tss(p_dict.get(x['Gene_ID'], []), x["enhancer_midpoint"]), axis=1)
            # # 过滤无法找到TSS的行
            # df1.dropna(subset=["promoter_data"], inplace=True)

            df1[["tss_chr", "tss", "tss_ID", "strand"]] = pd.DataFrame(df1["promoter_data"].tolist(), index=df1.index)
            df1.drop(columns=["enhancer_data", "promoter_data"], inplace=True)
            df1["enhancer_length"] = df1["enhancer_end"] - df1["enhancer_start"]
            df1['chr_size'] = df1['chr'].map(hg19_chromsize)
            df1["promoter_start"], df1["promoter_end"] = zip(*df1.apply(lambda x: compute_promoter_region(x["tss"], x["strand"], x['chr_size']), axis=1))

            # 计算距离
            df1["distance"] = df1.apply(lambda row: abs(row['enhancer_midpoint'] - row['tss']), axis=1)

            df_out = df1[df1["ccRE_type"] == "Enhancer-like"][["chr", "enhancer_start", "enhancer_end", "promoter_start", "promoter_end", "distance", "label", "cv", "strand", "ccRE_ID", "Gene_ID", "tss_ID", "ccRE_type", "enhancer_midpoint", "enhancer_length"]]
            
            # 筛选出正例集，即label为1的行
            positive_examples = df_out[df_out["label"] == 1].copy()

            # 计算当前文件的平均增强子长度
            current_avg_enhancer_length = positive_examples['enhancer_length'].mean()
            print(f'{data_file}: Sample size={len(positive_examples)}, Average Enhancer Length={current_avg_enhancer_length}')

            # 累加总样本量和总长度
            total_samples += len(positive_examples)
            total_enhancer_length += positive_examples['enhancer_length'].sum()

            # 删除不再需要的列
            positive_examples.drop(columns=["enhancer_length"], inplace=True)
            ep_output_file = os.path.join(output_dir, f'{os.path.splitext(data_file)[0]}_ep_distanceisWindows.txt')
            # positive_examples.to_csv(ep_output_file, sep="\t", index=None, header=True) #包含列名，以后就可以直接读取，不用再给列命名了
            process_and_save_samples(positive_examples, transcripts_df, ep_output_file)

    # 计算整体平均增强子长度
    overall_avg_enhancer_length = total_enhancer_length / total_samples if total_samples > 0 else 0
    print(f'Overall Average Enhancer Length: {overall_avg_enhancer_length}')

if __name__ == "__main__":
    main()


# '''统计距离'''
# import pandas as pd

# # 加载数据
# positive_examples = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio/NHEK.HiC-Benchmark.v3_ep_distanceisWindows.txt', sep='\t')

# # 设置 bin 的大小和范围
# bin_size = 2000
# max_distance = int(positive_examples['distance'].max()) + bin_size  # 确保包含最大值，且强制为整数

# # 生成 bins
# bins = range(0, max_distance, bin_size)

# # 划分 bins
# distance_bins = pd.cut(positive_examples['distance'], bins=bins, right=False)

# # 计算每个 bin 的样本数
# bin_counts = distance_bins.value_counts().sort_index()

# # 计算累积样本数
# cumulative_counts = bin_counts.cumsum()

# # 找到从 42000 开始累积样本数达到总样本的 80% 的最小 bin
# total_samples = positive_examples.shape[0]
# target_count = 0.8 * total_samples
# cumulative_counts_from_42000 = cumulative_counts[cumulative_counts.index >= pd.Interval(42000, 42000+bin_size, closed='left')]

# # 找到第一个使得累积样本数量从 42000 开始达到目标的 bin
# k_bin = cumulative_counts_from_42000[cumulative_counts_from_42000 >= target_count].index[0]

# # 打印结果
# print(f"Bin where cumulative count reaches at least 80% of total samples starting from 42000: {k_bin}")
# print(f"This corresponds to a distance from 42000 up to {k_bin.right}")


# '''筛选TSS'''
# # 加载 GTF 文件
# import pandas as pd
# import re
# from misc_utils import hg19_chromsize

# def compute_promoter_region(tss, strand, chr_size):
#     if strand == '+':
#         start = max(1, tss - 1500)
#         end = min(tss + 500, chr_size)
#     else:
#         start = max(1, tss - 500)
#         end = min(tss + 1500, chr_size)
#     return start, end

# def standardize_gene_id(gene_id):
#     # 正则表达式去除版本号和其他后缀
#     match = re.match(r"(ENSG\d+)", gene_id)
#     if match:
#         return match.group(1)
#     return None
# def standardize_transcript_id(transcript_id):
#     # 使用正则表达式去除版本号和其他后缀
#     match = re.match(r"(ENST\d+)", transcript_id)
#     if match:
#         return match.group(1)
#     return None

# # 加载 GTF 文件
# gtf_path = '/home/tjzhang03/zxj/deal_data/genomic_data/gencode.v43lift37.annotation.gtf'
# gtf_data = pd.read_csv(gtf_path, comment='#', sep='\t', header=None, names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# # 筛选出转录本
# transcripts = gtf_data[gtf_data['feature'] == 'transcript']

# # 解析属性字段
# def parse_attributes(attribute_str):
#     attributes = {}
#     for attr in attribute_str.split(';'):
#         if attr.strip():
#             key_value = attr.strip().split(' ')
#             if len(key_value) == 2:
#                 key, value = key_value
#                 attributes[key.strip()] = value.strip('"')
#     return attributes

# transcripts['attributes'] = transcripts['attribute'].apply(parse_attributes)
# transcripts = transcripts[(transcripts['attributes'].apply(lambda x: x.get('transcript_type') in ['protein_coding', 'lincRNA', 'antisense_RNA'])) &
#                           (transcripts['attributes'].apply(lambda x: int(x.get('level', 3)) <= 2))]

# # 提取 TSS 并计算启动子区域
# def compute_promoter(tss, strand):
#     return (max(1, tss - 1500), tss + 500) if strand == '+' else (tss - 500, max(1, tss + 1500))

# transcripts['TSS'] = transcripts.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)
# transcripts['chr_size'] = transcripts['chr'].map(hg19_chromsize)
# transcripts[['promoter_start', 'promoter_end']] = transcripts.apply(lambda x: compute_promoter_region(x['TSS'], x['strand'], x['chr_size']), axis=1, result_type='expand')

# # # 提取需要的信息
# # transcripts['Gene_ID'] = transcripts['attributes'].apply(lambda x: x.get('gene_id'))
# # transcripts['transcript_id'] = transcripts['attributes'].apply(lambda x: x.get('transcript_id'))
# # 标准化 gene_id 和 transcript_id
# transcripts['Gene_ID'] = transcripts['attributes'].apply(lambda x: standardize_gene_id(x.get('gene_id')))
# transcripts['tss_ID'] = transcripts['attributes'].apply(lambda x: standardize_transcript_id(x.get('transcript_id')))


# # 选择需要的列
# final_transcripts = transcripts[['chr', 'promoter_start', 'promoter_end', 'strand', 'TSS', 'Gene_ID', 'tss_ID']]

# # 显示结果
# print(final_transcripts.head())



# '''创造阴性样本'''
# import numpy as np

# # 加载数据
# positive_examples = pd.read_csv('/home/tjzhang03/zxj/deal_data/data_process/EG2EP/All-Pairs.Natural-Ratio/NHEK.HiC-Benchmark.v3_ep_distanceisWindows.txt', sep='\t')
# positive_examples['Gene_ID'] = positive_examples['Gene_ID'].apply(standardize_gene_id)
# positive_examples['tss_ID'] = positive_examples['tss_ID'].apply(standardize_transcript_id)
# # 筛选出距离范围内的样本
# pos_samples_df = positive_examples[(positive_examples['distance'] >= 42000) & (positive_examples['distance'] <= 500000)]
# pos_samples_df = pos_samples_df[['chr', 'enhancer_start', 'enhancer_end', 'promoter_start', 'promoter_end', 'distance', 'label', 'strand', 'ccRE_ID', 'Gene_ID', 'tss_ID', 'enhancer_midpoint']]
# # 预计算每个增强子的正样本数量
# sample_counts = pos_samples_df['ccRE_ID'].value_counts()

# # 获取增强子信息，避免后续重复的drop_duplicates操作
# positive_enhancers = pos_samples_df.drop_duplicates()
# # 使用set存储Gene_ID，用于快速查找和检查
# positive_genes = set(positive_enhancers['Gene_ID'])

# def generate_negative_samples(enhancer_row, transcripts_df, positive_genes, sample_counts):
#     # 筛选条件
#     mask = (
#         (transcripts_df['chr'] == enhancer_row['chr']) &
#         (~transcripts_df['Gene_ID'].isin(positive_genes)) &
#         (42000 <= np.abs(transcripts_df['TSS'] - enhancer_row['enhancer_midpoint'])) &
#         (np.abs(transcripts_df['TSS'] - enhancer_row['enhancer_midpoint']) <= 500000)
#     )
#     valid_transcripts = transcripts_df[mask]
#     num_pos_samples = sample_counts.get(enhancer_row['ccRE_ID'], 0)

#     # 检查可用的负样本数量并调整请求的样本数量
#     num_samples_to_return = min(len(valid_transcripts), num_pos_samples)
#     if num_samples_to_return < num_pos_samples:
#         print(f"Not enough negative samples for CCRE_ID {enhancer_row['ccRE_ID']}. Expected {num_pos_samples}, got {num_samples_to_return}.")
    
#     # 如果有足够的负样本候选，随机选择与正样本数量相同的负样本
#     if len(valid_transcripts) > 0:
#         samples = valid_transcripts.sample(n=num_samples_to_return)
#         samples['enhancer_start'] = enhancer_row['enhancer_start']
#         samples['enhancer_end'] = enhancer_row['enhancer_end']
#         samples['label'] = 0  # 设置负样本标签为0
#         samples['distance'] = np.abs(samples['TSS'] - enhancer_row['enhancer_midpoint'])
#         samples['ccRE_ID'] = enhancer_row['ccRE_ID']
#         return samples[['chr', 'enhancer_start', 'enhancer_end', 'promoter_start', 'promoter_end', 'distance', 'label', 'strand', 'ccRE_ID', 'Gene_ID', 'tss_ID']]
#     else:
#         return pd.DataFrame()
        
# # 应用筛选函数
# generated_samples = set()
# negative_samples = pd.DataFrame()  #避免对同一 CCRE_ID 多次生成负样本
# for index, enhancer in positive_enhancers.iterrows():
#     if enhancer['ccRE_ID'] not in generated_samples:
#         neg_samples = generate_negative_samples(enhancer, final_transcripts, positive_genes, sample_counts)
#         negative_samples = pd.concat([negative_samples, neg_samples], ignore_index=True)
#         generated_samples.add(enhancer['ccRE_ID'])

# # negative_samples = pd.DataFrame()
# # for index, enhancer in positive_enhancers.iterrows():
# #     neg_samples = generate_negative_samples(enhancer, final_transcripts, positive_genes)
# #     negative_samples = pd.concat([negative_samples, neg_samples], ignore_index=True)
# print(f"Number of Positive Samples: {pos_samples_df.shape[0]}")
# print(f"Number of Negative Samples: {negative_samples.shape[0]}")
# # 合并正负样本并按照ccRE_ID进行分组，确保每组内正样本紧跟其负样本
# all_samples = pd.concat([pos_samples_df, negative_samples], ignore_index=True)
# all_samples_sorted = all_samples.sort_values(by=['ccRE_ID', 'label'], ascending=[True, False])
# # 保存DataFrame到CSV
# all_samples_sorted = all_samples_sorted[['chr', 'enhancer_start', 'enhancer_end', 'promoter_start', 'promoter_end', 'distance', 'label', 'strand', 'ccRE_ID', 'Gene_ID', 'tss_ID']]
# output_file_path = '/home/tjzhang03/zxj/deal_data/data_process/create_data/NHEK.HiC-Benchmark.v3_ep_distanceisWindows.txt'
# all_samples_sorted.to_csv(output_file_path, index=False, sep='\t')