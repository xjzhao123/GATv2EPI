#!/bin/bash

# 定义文件数组
files=(
    "/home/tjzhang03/zxj/GATv2EPI/TransEPI/NHEK/NHEK.24_val.v3.tsv.gz"
    "/home/tjzhang03/zxj/GATv2EPI/TransEPI/NHEK/NHEK.24_train.v3.tsv.gz"
    "/home/tjzhang03/zxj/GATv2EPI/TransEPI/NHEK/NHEK.24_test.v3.tsv.gz"
)

# 循环处理每个文件
for file in "${files[@]}"; do
    # 获取文件路径和扩展名
    dir=$(dirname "$file")
    base=$(basename "$file" .gz) # 去掉.gz
    new_name="${base/24_/}" # 替换24_为空
    new_file="$dir/$new_name.tsv.gz" # 新的压缩文件名

    # 解压文件并重命名为去掉24_的新名称
    gunzip -c "$file" > "$dir/$new_name.tsv" # 直接解压并命名

    # 压缩成新的文件名
    gzip -c "$dir/$new_name.tsv" > "$new_file" # 压缩并输出到新文件
done



