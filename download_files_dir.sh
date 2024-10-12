#!/bin/bash

# 运行指令：
# chmod +x download_files.sh
# ./download_files.sh

# 定义CSV文件路径
CSV_FILE="info.csv"

# 定义本地下载目录
PDB_DOWNLOAD_DIR="pdb_download"
JSON_DOWNLOAD_DIR="json_download"

# 创建下载目录
mkdir -p $PDB_DOWNLOAD_DIR
mkdir -p $JSON_DOWNLOAD_DIR

# SSH密码 此过程不需要
# SSH_PASSWORD="-"

# 读取CSV文件并处理每一行
# 使用awk跳过标题行，并确保使用逗号作为分隔符
tail -n +2 "$CSV_FILE" | while IFS=',' read -r id batch_result_dir done_tag result_type target_pdb_name plddt ptm iptm pdockq

do
    # 定义远程文件路径
    PDB_FILE_PATH="$batch_result_dir/$target_pdb_name"
    
    # 构建JSON文件名
    JSON_FILE_NAME="${target_pdb_name/relaxed/scores}"
    JSON_FILE_PATH="$batch_result_dir/${JSON_FILE_NAME%.pdb}.json"

    # 下载PDB文件
    scp user@example.com:"$PDB_FILE_PATH" "$PDB_DOWNLOAD_DIR/"

    # 下载JSON文件
    scp user@example.com:"$JSON_FILE_PATH" "$JSON_DOWNLOAD_DIR/"

done
