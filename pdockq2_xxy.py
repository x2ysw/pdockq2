import numpy as np
import sys
import os
import argparse
import pickle
import pandas as pd
import itertools
from scipy.optimize import curve_fit
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

def retrieve_IFplddt(structure, chain1, chain2_lst, max_dist):
    chain_lst = list(chain1) + chain2_lst
    ifplddt = []
    contact_chain_lst = []
    for res1 in structure[0][chain1]:
        for chain2 in chain2_lst:
            count = 0
            for res2 in structure[0][chain2]:
                if res1.has_id('CA') and res2.has_id('CA'):
                    dis = abs(res1['CA'] - res2['CA'])
                    if dis <= max_dist:
                        ifplddt.append(res1['CA'].get_bfactor())
                        count += 1
                elif res1.has_id('CB') and res2.has_id('CB'):
                    dis = abs(res1['CB'] - res2['CB'])
                    if dis <= max_dist:
                        ifplddt.append(res1['CB'].get_bfactor())
                        count += 1
            if count > 0:
                contact_chain_lst.append(chain2)
    contact_chain_lst = sorted(list(set(contact_chain_lst)))
    if len(ifplddt) > 0:
        IF_plddt_avg = np.mean(ifplddt)
    else:
        IF_plddt_avg = 0
    return IF_plddt_avg, contact_chain_lst

def retrieve_IFPAEinter(structure, paeMat, contact_lst, max_dist):
    paeMat = np.array(paeMat)
    chain_lst = [x.id for x in structure[0]]
    seqlen = [len(x) for x in structure[0]]
    ifpae_avg = []
    d = 10
    for ch1_idx in range(len(chain_lst)):
        idx = chain_lst.index(chain_lst[ch1_idx])
        ch1_sta = sum(seqlen[:idx])
        ch1_end = ch1_sta + seqlen[idx]
        ifpae_col = []
        for contact_ch in contact_lst[ch1_idx]:
            index = chain_lst.index(contact_ch)
            ch_sta = sum(seqlen[:index])
            ch_end = ch_sta + seqlen[index]
            remain_paeMatrix = paeMat[ch1_sta:ch1_end, ch_sta:ch_end]
            mat_x = -1
            for res1 in structure[0][chain_lst[ch1_idx]]:
                mat_x += 1
                mat_y = -1
                for res2 in structure[0][contact_ch]:
                    mat_y += 1
                    if res1.has_id('CA') and res2.has_id('CA'):
                        dis = abs(res1['CA'] - res2['CA'])
                        if dis <= max_dist:
                            ifpae_col.append(remain_paeMatrix[mat_x, mat_y])
        if not ifpae_col:
            ifpae_avg.append(0)
        else:
            norm_if_interpae = np.mean(1 / (1 + (np.array(ifpae_col) / d) ** 2))
            ifpae_avg.append(norm_if_interpae)
    return ifpae_avg

def calc_pmidockq(ifpae_norm, ifplddt):
    df = pd.DataFrame()
    df['ifpae_norm'] = ifpae_norm
    df['ifplddt'] = ifplddt
    df['prot'] = df.ifpae_norm * df.ifplddt
    fitpopt = [1.31034849e+00, 8.47326239e+01, 7.47157696e-02, 5.01886443e-03]  # from original fit function
    df['pmidockq'] = sigmoid(df.prot.values, *fitpopt)
    return df

def sigmoid(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y

def process_files(pdb_file, pkl_file, dist):
    # 解析pdb文件
    pdbp = PDBParser(QUIET=True)
    structure = pdbp.get_structure('', pdb_file)
    # 获取链信息
    chains = []
    for chain in structure[0]:
        chains.append(chain.id)
    # 计算界面plDDT值
    remain_contact_lst = []
    # retrieve interface plDDT at chain-level
    plddt_lst = []
    for idx in range(len(chains)):
        chain2_lst = list(set(chains) - set(chains[idx]))
        IF_plddt, contact_lst = retrieve_IFplddt(structure, chains[idx], chain2_lst, dist)
        plddt_lst.append(IF_plddt)
        remain_contact_lst.append(contact_lst)
    # 读取pickle文件并计算界面PAE值
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    avgif_pae = retrieve_IFPAEinter(structure, data['predicted_aligned_error'], remain_contact_lst, dist)
    # calculate pmiDockQ
    res = calc_pmidockq(avgif_pae, plddt_lst)
    # 返回结果
    return chains, res['pmidockq'].tolist()

def main():
    parser = argparse.ArgumentParser(description='Calculate chain_level pDockQ_i in batch.')
    
    parser.add_argument('--pdb_dir', nargs=1, type=str, required=True, help='Directory containing pdb files.')
    parser.add_argument('--pkl_dir', nargs=1, type=str, required=True, help='Directory containing pkl files.')
    parser.add_argument('--input_csv', nargs=1, type=str, required=True, help='Input CSV file with pdb and pkl filenames.')
    parser.add_argument('--output_csv', nargs=1, type=str, required=True, help='Output CSV file name.')
    parser.add_argument("-dist", help="Maximum distance of a contact", nargs='?', type=int, default=8)
    args = parser.parse_args()
    print(args)

    # 读取CSV文件
    df_files = pd.read_csv(args.input_csv[0])

    # 初始化结果列表
    results = []

    # 遍历每一行，处理文件
    for index, row in df_files.iterrows():
        pdb_filename = row[0]
        pkl_filename = row[1]
        pdb_file = os.path.join(args.pdb_dir[0], pdb_filename)
        pkl_file = os.path.join(args.pkl_dir[0], pkl_filename)

        print(f'Processing PDB file: {pdb_filename}, PKL file: {pkl_filename}')

        # 检查文件是否存在
        if not os.path.exists(pdb_file):
            print(f'PDB file {pdb_file} does not exist.')
            continue
        if not os.path.exists(pkl_file):
            print(f'PKL file {pkl_file} does not exist.')
            continue

        try:
            chains, scores = process_files(pdb_file, pkl_file, args.dist)
            # 打印每个链的得分
            for ch, score in zip(chains, scores):
                print(f'Chain {ch}: Score {score}')
            # 将结果添加到列表
            result_row = [pdb_filename] + scores
            results.append(result_row)
        except Exception as e:
            print(f'Error processing files {pdb_filename} and {pkl_filename}: {e}')
            continue

    # 确定最大链数量（用于CSV列数）
    max_chains = max(len(r) - 1 for r in results)

    # 准备CSV列名
    columns = ['PDB_File']
    for i in range(1, max_chains + 1):
        columns.append(f'Chain_{i}_Score')

    # 创建DataFrame并保存为CSV
    df_results = pd.DataFrame(results, columns=columns)
    df_results.to_csv(args.output_csv[0], index=False)
    print(f'Results saved to {args.output_csv[0]}')

if __name__ == '__main__':
    main()
