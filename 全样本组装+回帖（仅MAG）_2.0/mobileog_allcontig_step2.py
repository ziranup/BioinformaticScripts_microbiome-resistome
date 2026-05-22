#!/usr/bin/env python3
# 脚本名称：MGEs_step2.py
# 功能：基于 DIAMOND 注释结果和原始 reads，计算每个样本中每个 MGE 的丰度（CPM、TPM、RPKM）
# 输入：
#   - DIAMOND 注释结果：./MGEs_mobileOG/mobileOG/result/mobileog_hits.tsv
#   - 非冗余基因集：./D_7_nonredundantGenes/nonredundant_genes.fasta
#   - clean reads目录：./D_2_HRreads/cleanreads/
# 输出：./MGEs_mobileOG/samtool/result/ 中的丰度矩阵文件

import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess
import argparse
from pathlib import Path
import re

def create_directories():
    """创建必要的输出目录"""
    os.makedirs("./MGEs_mobileOG/samtool/result", exist_ok=True)
    os.makedirs("./MGEs_mobileOG/samtool/tmp", exist_ok=True)

def parse_diamond_results(diamond_file):
    """
    解析 DIAMOND blastp 结果文件
    返回 DataFrame 包含 gene_id -> mobileOG_id 的映射
    """
    print(f"解析 DIAMOND 结果文件: {diamond_file}")
    
    # 读取 DIAMOND 结果
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'stitle']
    
    # 添加 low_memory=False 避免类型推断问题
    df = pd.read_csv(diamond_file, sep='\t', names=columns, comment='#', low_memory=False)
    
    # 强制转换数值列的数据类型
    numeric_columns = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 
                      'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    
    for col in numeric_columns:
        if col in df.columns:
            # 使用 errors='coerce' 将无法转换的值设为 NaN
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # 删除包含 NaN 的行（这些是转换失败的行）
    original_len = len(df)
    df = df.dropna(subset=numeric_columns)
    if len(df) < original_len:
        print(f"警告: 删除了 {original_len - len(df)} 行包含无效数值的记录")
    
    # 过滤低质量比对
    df = df[(df['pident'] >= 80) & (df['evalue'] <= 1e-6)]
    
    # 提取 mobileOG ID 和功能信息
    df['mobileOG_id'] = df['sseqid']
    df['function'] = df['stitle'].str.split(' ', n=1).str[1].fillna('')
    
    # 只保留需要的列
    annotation_df = df[['qseqid', 'mobileOG_id', 'function']].drop_duplicates()
    
    print(f"共注释到 {len(annotation_df)} 个非冗余基因")
    return annotation_df

def get_gene_lengths(fasta_file):
    """
    从 FASTA 文件获取基因长度
    """
    print(f"读取基因长度信息: {fasta_file}")
    gene_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_lengths[record.id] = len(record.seq)
    return gene_lengths

def get_sample_list(cleanreads_dir):
    """
    获取样本列表（从 clean reads 目录）
    """
    r1_files = list(Path(cleanreads_dir).glob("*hrm.1.fastq"))
    samples = []
    for f in r1_files:
        sample_name = f.name.replace('.hrm.1.fastq', '')
        samples.append(sample_name)
    return sorted(samples)

def count_reads_in_file(file_path):
    """
    计算 FASTQ 文件中的 reads 数量
    """
    result = subprocess.run(['wc', '-l', str(file_path)], capture_output=True, text=True)
    lines = int(result.stdout.split()[0])
    return lines // 4

def map_reads_to_genes(sample_name, cleanreads_dir, nonredundant_fasta, tmp_dir):
    """
    使用 bowtie2 将 reads 映射到非冗余基因集，并统计每个基因的 reads 数
    """
    print(f"处理样本: {sample_name}")
    
    # 构建 bowtie2 索引（如果不存在）
    index_prefix = os.path.join(tmp_dir, "nonredundant_genes_index")
    if not os.path.exists(f"{index_prefix}.1.bt2"):
        print("构建 bowtie2 索引...")
        # 修改这里：替换 capture_output=True
        subprocess.run([
            'bowtie2-build', 
            nonredundant_fasta, 
            index_prefix
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # 获取 reads 文件路径
    r1_file = os.path.join(cleanreads_dir, f"{sample_name}.hrm.1.fastq")
    r2_file = os.path.join(cleanreads_dir, f"{sample_name}.hrm.2.fastq")
    
    # 创建 SAM 文件路径
    sam_file = os.path.join(tmp_dir, f"{sample_name}.sam")
    
    # 运行 bowtie2 映射
    bowtie2_cmd = [
        'bowtie2',
        '-x', index_prefix,
        '-1', r1_file,
        '-2', r2_file,
        '-S', sam_file,
        '--very-sensitive',
        '--no-unal',
        '-p', '8'
    ]
    subprocess.run(bowtie2_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # 统计每个基因的 reads 数
    gene_counts = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                gene_id = fields[2]
                if gene_id != '*':
                    gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1
    
    # 清理 SAM 文件
    os.remove(sam_file)
    
    return gene_counts


def calculate_abundances(annotation_df, gene_lengths, sample_gene_counts, samples):
    """
    计算 CPM、TPM、RPKM 丰度
    """
    print("计算丰度指标...")
    
    # 创建 MGE 丰度矩阵
    mobileOG_ids = sorted(annotation_df['mobileOG_id'].unique())
    abundance_matrix_cpm = pd.DataFrame(0, index=mobileOG_ids, columns=samples, dtype=float)
    abundance_matrix_tpm = pd.DataFrame(0, index=mobileOG_ids, columns=samples, dtype=float)
    abundance_matrix_rpkm = pd.DataFrame(0, index=mobileOG_ids, columns=samples, dtype=float)
    
    # 为每个样本计算丰度
    for sample in samples:
        print(f"  处理样本 {sample}...")
        gene_counts = sample_gene_counts[sample]
        
        # 将基因 counts 映射到 MGE
        mge_counts = {}
        for _, row in annotation_df.iterrows():
            gene_id = row['qseqid']
            mobileOG_id = row['mobileOG_id']
            if gene_id in gene_counts:
                mge_counts[mobileOG_id] = mge_counts.get(mobileOG_id, 0) + gene_counts[gene_id]
        
        # 计算总 mapped reads
        total_mapped_reads = sum(gene_counts.values())
        if total_mapped_reads == 0:
            continue
        
        # 计算 CPM (Counts Per Million)
        for mobileOG_id, count in mge_counts.items():
            cpm = (count / total_mapped_reads) * 1e6
            abundance_matrix_cpm.loc[mobileOG_id, sample] = cpm
        
        # 计算 RPKM/FPKM (需要基因长度信息)
        # 对于每个 MGE，我们需要其对应的基因长度
        mge_lengths = {}
        for _, row in annotation_df.iterrows():
            gene_id = row['qseqid']
            mobileOG_id = row['mobileOG_id']
            if gene_id in gene_lengths and mobileOG_id not in mge_lengths:
                mge_lengths[mobileOG_id] = gene_lengths[gene_id]
        
        # 计算 RPKM
        for mobileOG_id, count in mge_counts.items():
            if mobileOG_id in mge_lengths:
                gene_length_kb = mge_lengths[mobileOG_id] / 1000.0
                rpkm = (count / (gene_length_kb * total_mapped_reads)) * 1e9
                abundance_matrix_rpkm.loc[mobileOG_id, sample] = rpkm
        
        # 计算 TPM (Transcripts Per Million)
        # TPM = (RPKM / sum(RPKM)) * 1e6
        rpkm_values = abundance_matrix_rpkm[sample]
        rpkm_sum = rpkm_values.sum()
        if rpkm_sum > 0:
            tpm_values = (rpkm_values / rpkm_sum) * 1e6
            abundance_matrix_tpm[sample] = tpm_values
    
    return abundance_matrix_cpm, abundance_matrix_tpm, abundance_matrix_rpkm

def main():
    create_directories()
    
    # 文件路径定义
    DIAMOND_FILE = "./MGEs_mobileOG/mobileOG/result/mobileog_hits.tsv"
    NONREDUNDANT_FASTA = "./D_7_nonredundantGenes/nonredundant_genes.fasta"
    CLEANREADS_DIR = "./D_2_HRreads/cleanreads"
    TMP_DIR = "./MGEs_mobileOG/samtool/tmp"
    
    # 步骤1: 解析 DIAMOND 结果
    annotation_df = parse_diamond_results(DIAMOND_FILE)
    
    # 步骤2: 获取基因长度
    gene_lengths = get_gene_lengths(NONREDUNDANT_FASTA)
    
    # 步骤3: 获取样本列表
    samples = get_sample_list(CLEANREADS_DIR)
    print(f"发现 {len(samples)} 个样本: {samples}")
    
    # 步骤4: 为每个样本映射 reads 并统计基因 counts
    sample_gene_counts = {}
    for sample in samples:
        gene_counts = map_reads_to_genes(sample, CLEANREADS_DIR, NONREDUNDANT_FASTA, TMP_DIR)
        sample_gene_counts[sample] = gene_counts
    
    # 步骤5: 计算丰度
    cpm_matrix, tpm_matrix, rpkm_matrix = calculate_abundances(
        annotation_df, gene_lengths, sample_gene_counts, samples
    )
    
    # 步骤6: 保存结果
    cpm_matrix.to_csv("./MGEs_mobileOG/samtool/result/mobileOG_abundance_CPM.tsv", sep='\t')
    tpm_matrix.to_csv("./MGEs_mobileOG/samtool/result/mobileOG_abundance_TPM.tsv", sep='\t')
    rpkm_matrix.to_csv("./MGEs_mobileOG/samtool/result/mobileOG_abundance_RPKM.tsv", sep='\t')
    
    # 保存注释信息
    annotation_df.to_csv("./MGEs_mobileOG/samtool/result/mobileOG_annotation.tsv", sep='\t', index=False)
    
    print("✅ 丰度计算完成！")
    print("结果文件:")
    print("  - CPM: ./MGEs_mobileOG/samtool/result/mobileOG_abundance_CPM.tsv")
    print("  - TPM: ./MGEs_mobileOG/samtool/result/mobileOG_abundance_TPM.tsv")
    print("  - RPKM: ./MGEs_mobileOG/samtool/result/mobileOG_abundance_RPKM.tsv")
    print("  - 注释: ./MGEs_mobileOG/samtool/result/mobileOG_annotation.tsv")

if __name__ == "__main__":
    main()