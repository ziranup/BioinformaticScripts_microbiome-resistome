#!/bin/bash
# 脚本名称：6_prodigal_orf.sh
# 功能：从组装 contigs 预测 ORFs（CDS），用于后续功能注释
# 前置数据准备：
#   - contigs.fa（来自 MEGAHIT 或 SPAdes）
# 输出：
#   - output_dir/genes.gff  # 包含ORF的 基因组坐标信息（起始位点、终止位点、链方向等）
#   - output_dir/genes.fna  # 包含ORF的核苷酸序列（即 DNA 序列）
#   - output_dir/proteins.faa  # 包含ORF的蛋白质序列（即氨基酸序列）

set -e

output_dir="./D_6_ORFs"
mkdir -p "${output_dir}"

prodigal -i ./D_3_contigs/mg_6.contigs.fa \
         -o "${output_dir}/genes.gff" \
         -d "${output_dir}/genes.fna" \
         -a "${output_dir}/proteins.faa" \
         -p meta \
         -f gff

echo "✅ Prodigal ORF 预测完成"