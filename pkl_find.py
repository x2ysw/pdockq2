import os
import json
import pickle

# 定义 JSON 文件所在目录
json_dir = '......./json_download'
# 定义 PKL 文件保存目录
pkl_dir = '......../pkl_turn'

# 确保 PKL 文件保存目录存在
os.makedirs(pkl_dir, exist_ok=True)

# 遍历目录中的所有文件
for filename in os.listdir(json_dir):
    if filename.endswith('.json'):
        # 构建完整的 JSON 文件路径
        json_file_path = os.path.join(json_dir, filename)
        
        # 读取 JSON 文件
        with open(json_file_path, 'r') as json_file:
            data = json.load(json_file)
        
        # 提取 PAE 矩阵
        pae_matrix = data.get('pae', None)
        
        # 构建要保存到 .pkl 文件中的数据结构
        pkl_data = {
            'predicted_aligned_error': pae_matrix
        }
        
        # 构建 PKL 文件路径
        pkl_file_path = os.path.join(pkl_dir, filename.replace('.json', '.pkl'))
        
        # 保存到 .pkl 文件
        with open(pkl_file_path, 'wb') as pkl_file:
            pickle.dump(pkl_data, pkl_file)
        
        print(f"Data has been successfully converted and saved to '{pkl_file_path}'.")

print("All JSON files have been successfully converted to PKL files.")
