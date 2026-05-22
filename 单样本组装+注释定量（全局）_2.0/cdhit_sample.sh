#!/bin/bash
# 脚本名称：10_cdhit_nonredundant_per_sample.sh
# 功能：为每个样本单独构建非冗余基因集（用于后续分析）
# 前置数据准备：
#   - ./D_a2_SampleORFs/{Sample}_genes.fna（来自 Prodigal 的单样本 ORF 核苷酸序列）
# 输出：
#   - ./D_7_nonredundantGenes_perSample/{Sample}_nonredundant_genes.fna

# cdhit参数说明（cd-hit-est 用于核苷酸序列）
# -i：fasta格式的输入序列文件
# -o：输出文件的文件名前缀
# -c：序列相似度identity阈值，默认为0.9。此处设为0.95以去除高度相似的冗余ORF
# -G：全局比对（默认1），适用于完整或接近完整的基因
# -M：内存限制(MB)，设为0则无限制
# -T：使用的CPU线程数，设为12以加速单个样本处理（可根据服务器负载调整）
# -d：聚类信息中序列名长度，设为0保留完整ID
# -n：word长度，对于-c 0.95，推荐-n 10（但cd-hit-est在-c>0.9时自动优化，可省略）
# -g 1：开启精确模式，确保最佳聚类（推荐用于构建参考基因集）

set -e

# 1. 创建输出目录
output_dir="./D_a3_SampleNonredundantGenes"
mkdir -p "${output_dir}"

# 2. 定义输入目录
input_orf_dir="./D_a2_SampleORFs"

# 3. 获取所有样本的唯一名称列表
# 通过 *_genes.fna 文件提取样本名 (如 CK_1, T1_1 等)
samples=$(ls ${input_orf_dir}/*_genes.fna | xargs -n1 basename | sed 's/_genes\.fna//')

# 4. 遍历每个样本，进行独立的非冗余化
for sample in $samples
do
    echo "🚀 正在处理样本: ${sample} ..."

    input_fna="${input_orf_dir}/${sample}_genes.fna"
    output_prefix="${output_dir}/${sample}_nonredundant_genes"

    # 运行 cd-hit-est 对单个样本的 genes.fna 进行去冗余
    cd-hit-est \
        -i "${input_fna}" \
        -o "${output_prefix}.fna" \
        -M 0 \
        -T 72 \
        -d 0 

    echo "✅ 样本 ${sample} 的非冗余基因集已生成"
done

echo "🎉 所有样本的非冗余基因集构建完成！"
echo "结果保存在: ./${output_dir}/"