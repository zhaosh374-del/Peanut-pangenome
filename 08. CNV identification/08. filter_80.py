#!/usr/bin/env python3
import os
from Bio import SeqIO
from collections import defaultdict

# 定义路径
protein_fasta = "/ref/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa"
input_dir = "/path/blastp_chr_same_uniq"
output_dir = "/path/blastp_chr_same_uniq_80_80"

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

def load_protein_lengths(fasta_file):
    """读取蛋白质FASTA文件并返回{protein_id: length}字典"""
    protein_lengths = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            protein_id = record.id  # 保留完整ID
            protein_lengths[protein_id] = len(record.seq)
    return protein_lengths

def process_blastp_file(input_path, output_path, protein_lengths):
    """处理单个BLASTP结果文件"""
    # 存储所有记录
    all_records = []
    # 存储Tifrunner ID到记录的映射
    tifrunner_to_records = defaultdict(list)
    
    with open(input_path, "r") as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 12:
                continue
                
            query_id = parts[0]
            subject_id = parts[1]
            identity = float(parts[2])
            alignment_length = int(parts[3])
            
            record = {
                'line': line,
                'query_id': query_id,
                'subject_id': subject_id,
                'identity': identity,
                'alignment_length': alignment_length,
                'protein_length': protein_lengths.get(subject_id, 0)
            }
            
            all_records.append(record)
            tifrunner_to_records[subject_id].append(record)
    
    # 确定哪些Tifrunner ID有多个MP ID映射
    multi_mapped = {tid for tid, recs in tifrunner_to_records.items() if len(recs) > 1}
    
    # 筛选记录
    filtered_records = []
    for record in all_records:
        subject_id = record['subject_id']
        
        # 如果是1:1映射，直接保留
        if subject_id not in multi_mapped:
            filtered_records.append(record)
            continue
            
        # 如果是多对一映射，应用过滤条件
        identity = record['identity']
        alignment_length = record['alignment_length']
        protein_length = record['protein_length']
        
        # 检查条件1：identity > 80
        if identity <= 80:
            continue
            
        # 检查条件2：alignment_length > 80%的蛋白质长度
        if protein_length > 0 and alignment_length >= 0.8 * protein_length:
            filtered_records.append(record)
    
    # 写入结果
    with open(output_path, "w") as outfile:
        for record in filtered_records:
            outfile.write(record['line'])

def main():
    # 第一步：加载蛋白质长度信息
    print("Loading protein lengths...")
    protein_lengths = load_protein_lengths(protein_fasta)
    print(f"Loaded lengths for {len(protein_lengths)} proteins")
    
    # 第二步：处理每个BLASTP结果文件
    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            
            print(f"Processing {filename}...")
            process_blastp_file(input_path, output_path, protein_lengths)
            print(f"Finished processing {filename}")
    
    print("All files processed successfully!")

if __name__ == "__main__":
    main()