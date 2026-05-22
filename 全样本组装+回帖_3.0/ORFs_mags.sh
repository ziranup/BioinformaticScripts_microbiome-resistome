#!/bin/bash

# 设置输入和输出路径
INPUT_DIR="./D_b5_MAGs/filtered_mags"
OUTPUT_DIR="./D_b6_ORFs_MAGs"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 定义单个样本的处理函数
# 此函数接受一个 fasta 文件路径作为参数
run_prodigal_on_mag() {
    local contig_file="$1"
    # 提取文件名（不含路径和扩展名），例如 bin.10.fa -> bin.10
    local sample_name=$(basename "$contig_file" .fa)
    
    # 运行 Prodigal
    # -i: 输入文件
    # -o: GFF 输出 (坐标信息)
    # -d: 核苷酸序列输出 (CDS)
    # -a: 氨基酸序列输出 (蛋白)
    # -p meta: 宏基因组模式 (适用于环境样本或单细胞)
    # -f gff: 输出格式
    # -q: 静默模式
    prodigal -i "$contig_file" \
             -o "${OUTPUT_DIR}/${sample_name}_genes.gff" \
             -d "${OUTPUT_DIR}/${sample_name}_genes.fna" \
             -a "${OUTPUT_DIR}/${sample_name}_proteins.faa" \
             -p meta \
             -f gff \
             -q
}

# 导出函数，以便 parallel 调用
export -f run_prodigal_on_mag
export OUTPUT_DIR
export INPUT_DIR


# 获取所有 .fa 文件列表
FA_FILES=("$INPUT_DIR"/*.fa)

# 使用 GNU parallel 并行处理所有 MAGs
# 假设你有 72 个核心可用，这里使用 -j 36 控制并发数（避免 I/O 瓶颈）
parallel -j 36 run_prodigal_on_mag {} ::: "${FA_FILES[@]}"

echo "所有 MAGs 的 ORF 预测已完成，结果保存在 $OUTPUT_DIR"