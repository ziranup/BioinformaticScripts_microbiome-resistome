#!/bin/bash
# 脚本功能：批量对单样本组装的 contigs 进行 ORF 预测
# 输入目录：./D_a1_SampleContigs/SContigs
# 输出目录：./D_a2_SampleORFs

# 注意：fna,ffn,faa都属于fasta(.fa)格式：
# fna （fasta nucleic acid file）所有核酸序列信息
# ffn （fasta nucleotide coding regions file）所有基因的核酸序列信息
# faa （fasta Amino Acid file） 即所有基因对应的蛋白质序列信息

set -e

# 1. 定义路径
INPUT_DIR="./D_a1_SampleContigs/SContigs"
OUTPUT_DIR="./D_a2_SampleORFs"
mkdir -p "${OUTPUT_DIR}"

# 2. 遍历输入目录下的所有 contig 文件
# 假设文件名格式为 SampleName.contigs.fa
for contig_file in ${INPUT_DIR}/*.contigs.fa
do
    # 提取样本名 (例如从 CK_1.contigs.fa 提取 CK_1)
    sample_name=$(basename "$contig_file" .contigs.fa)
    
    echo "🚀 正在处理样本: ${sample_name} ..."

    # 3. 运行 Prodigal
    # -i: 输入文件
    # -o: GFF 输出 (坐标信息)
    # -d: 核苷酸序列输出 (CDS)
    # -a: 氨基酸序列输出 (蛋白)
    # -p meta: 宏基因组模式 (适用于环境样本或单细胞)
    # -f gff: 输出格式
    prodigal -i "$contig_file" \
             -o "${OUTPUT_DIR}/${sample_name}_genes.gff" \
             -d "${OUTPUT_DIR}/${sample_name}_genes.fna" \
             -a "${OUTPUT_DIR}/${sample_name}_proteins.faa" \
             -p meta \
             -f gff \
             -q
    
    echo "✅ ${sample_name} 预测完成"
done

echo "🎉 所有样本 ORF 预测已全部完成！"