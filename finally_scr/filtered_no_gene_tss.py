'''删除非基因'''
import pandas as pd
import re

def standardize_gene_id(gene_id):
    '''标准化基因ID ，有很大后缀，包括版本号之类的'''
    match = re.match(r"(ENSG\d+)", gene_id)
    if match:
        return match.group(1)
    return None

def load_gene_ids_from_gtf(gtf_file):
    # 读取GTF文件，跳过注释行，没有头部信息
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, usecols=[2, 8],
                      names=['type', 'attributes'])
    # 仅选择基因行
    gtf_genes = gtf[gtf['type'] == 'gene']
    # 提取基因ID
    gtf_genes['gene_id'] = gtf_genes['attributes'].str.extract(r'gene_id "(ENSG\d+)')
    # 标准化基因ID
    gtf_genes['gene_id'] = gtf_genes['gene_id'].apply(standardize_gene_id)
    return set(gtf_genes['gene_id'].dropna())

def filter_bed_by_genes(bed_file, gene_ids, output_file, non_gene_ids_file):
    # 读取BED文件
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed['standard_gene_id'] = bed.iloc[:, 6].apply(standardize_gene_id)
    # 计算初始行数
    initial_count = len(bed)
    # 过滤基因ID
    filtered_bed = bed[bed['standard_gene_id'].isin(gene_ids)]
    # 保存过滤后的文件
    filtered_bed.to_csv(output_file, sep='\t', index=False, header=False)
    # 计算并保存非基因ID
    non_gene_ids = bed[~bed['standard_gene_id'].isin(gene_ids)]['standard_gene_id'].dropna().unique()
    with open(non_gene_ids_file, 'w') as f:
        for id in non_gene_ids:
            f.write(id + '\n')
    # 计算删除的行数
    removed_count = initial_count - len(filtered_bed)
    return len(filtered_bed), removed_count, len(non_gene_ids)

# 文件路径
gtf_file = '/home/tjzhang03/zxj/deal_data/genomic_data/gencode.v43lift37.annotation.gtf'
bed_file = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/GENCODEv19-TSSs.bed'
output_file = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/GENCODEv19-TSSs_filtered.bed'
non_gene_ids_file = '/home/tjzhang03/zxj/deal_data/data_input/BENGI/non_gene_ids.txt'

# 加载有效的基因ID
gene_ids = load_gene_ids_from_gtf(gtf_file)
# 过滤.bed文件并记录非基因ID
kept_count, removed_count, non_gene_count = filter_bed_by_genes(bed_file, gene_ids, output_file, non_gene_ids_file)

print(f"Total lines kept: {kept_count}")
print(f"Total lines removed: {removed_count}")
print(f"Total non-gene IDs recorded: {non_gene_count}")