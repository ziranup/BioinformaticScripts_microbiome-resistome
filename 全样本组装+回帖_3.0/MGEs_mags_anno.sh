#!/bin/bash

# 使用 mobileOG 数据库注释 MAGs 中的移动遗传元件（MGEs）的脚本（2026年5月）
# 作者：YZR
# 说明：基于 DIAMOND blastp 比对 MAGs 预测的蛋白序列到 mobileOG 数据库

# 1. 设置 mobileOG 数据库路径（Diamond 格式）
MOBILEOG_DB="/database/work/zryan/biodb/mobileOG-db_1.6/mobileOG-db_beatrix-1.6.All.dmnd"

# 2. 输入：MAGs 预测的蛋白质文件目录（faa 文件）
PROTEINS_DIR="./D_b6_ORFs_MAGs"

# 3. 输出目录：存放每个 bin 的 MGE 注释结果
OUTPUT_DIR="./MGEs_mobileOG_MAGs"
mkdir -p "$OUTPUT_DIR"

# 4. 临时目录：用于 DIAMOND 比对过程中的临时文件（避免 /tmp 空间不足）
TMPDIR="./tmp_diamond_mobileog"
mkdir -p "$TMPDIR"

# 5. 定义处理单个样本的函数（关键：封装含 . 的命令，避免 parallel 报错）
annotate_mge() {
    local protein_file="$1"
    
    # 提取样本名，例如 bin.100_proteins.faa → bin.100
    # 注意：使用 _proteins.faa 作为后缀进行截断
    local sample_name=$(basename "$protein_file" _proteins.faa)
    
    # 运行 DIAMOND blastp 比对到 mobileOG 数据库
    diamond blastp \
        --db "$MOBILEOG_DB" \
        --query "$protein_file" \
        --out "$OUTPUT_DIR/${sample_name}_mobileOG.tsv" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
        --threads 2 \
        --evalue 1e-7 \
        --max-target-seqs 1 \
        --sensitive \
        --tmpdir "$TMPDIR" \
        --quiet \
        --header simple
}

# 6. 导出函数和所需变量，使 GNU Parallel 能在子 shell 中调用
export -f annotate_mge
export MOBILEOG_DB OUTPUT_DIR TMPDIR

# 7. 获取所有蛋白质文件列表（匹配 *_proteins.faa）
# 使用数组存储，避免空格或特殊字符问题
protein_files=("$PROTEINS_DIR"/*_proteins.faa)

# 检查是否找到文件
if [ ${#protein_files[@]} -eq 0 ] || [ ! -e "${protein_files[0]}" ]; then
    echo "错误：未在 $PROTEINS_DIR 中找到任何 *_proteins.faa 文件！"
    exit 1
fi

# 8. 并行调用函数处理所有样本
# -j 36：使用 36 个并行任务（可根据服务器负载调整，建议不超过总核心数的 1/2）
parallel -j 36 annotate_mge ::: "${protein_files[@]}"

echo "mobileOG 移动遗传元件（MGEs）注释完成！结果位于 $OUTPUT_DIR"