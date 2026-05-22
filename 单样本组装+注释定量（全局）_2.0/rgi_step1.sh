#!/bin/bash

# 脚本名称: run_rgi_annotation.sh
# 功能描述: 基于单样本组装的contigs，使用RGI+CARD注释ARGs（仅注释步骤，不含丰度计算）

# 说明：coverm的参考是contig而不是ORF.因此，rgi main不要使用-t protein模式，这样会导致后面的定量很难且不准确。这里使用-t contig会重复计算一次prodigal生成orf，但是可以得到更详细的注释结果并且方便定量。。

# 设置进程数（每个样本分配12线程，并行6个样本，总72线程）
THREADS_PER_SAMPLE=12
PARALLEL_JOBS=6

# 设置路径
# contig用于提供样本级别长序列
CONTIGS_DIR="./D_a1_SampleContigs/SContigs"

# 数据库路径
RGI_DB_PATH="/database/work/zryan/biodb/CARD/card.json"

# 输出目录
RGI_RESULT_DIR="./ARGs_CARD/rgi/result"
RGI_TMP_DIR="./ARGs_CARD/rgi/tmp"

# 创建输出目录
mkdir -p "$RGI_RESULT_DIR" "$RGI_TMP_DIR"

# === 检查 rgi 是否可用 ===
if ! command -v rgi &> /dev/null; then
    echo "错误: 'rgi' 命令未找到。"
    echo "请确保已激活包含 rgi (from CARD) 的 conda 环境。"
    echo "安装建议：mamba create -n rgi -c conda-forge -c bioconda rgi diamond"
    exit 1
fi

# 导入RGI数据库（只需一次）
echo "正在加载CARD数据库到RGI..."
rgi load --card_json "$RGI_DB_PATH" --local

# 提取所有样本名（从contigs文件名推断）
# 示例：CK_1.contigs.fa → CK_1
SAMPLES=$(ls "$CONTIGS_DIR"/*.contigs.fa | xargs -n1 basename | sed 's/\.contigs\.fa$//')

# 定义处理单个样本的函数
process_sample() {
    local sample=$1
    local threads=$2

    # 路径定义
    local contig_file="$CONTIGS_DIR/${sample}.contigs.fa"

    # RGI输出
    local rgi_out_prefix="$RGI_RESULT_DIR/${sample}_arg"
    local rgi_txt="${rgi_out_prefix}.txt"

    # 步骤: RGI注释（contig模式）
    if [ ! -f "$rgi_txt" ] || [ ! -s "$rgi_txt" ]; then
        echo "正在处理样本: $sample"
        rgi main \
            -i "$contig_file" \
            -o "$rgi_out_prefix" \
            -t contig \
            -a DIAMOND \
            -n "$threads" \
            --include_nudge \
            --local \
            --debug
        echo "✅ 样本 $sample RGI注释完成"
    else
        echo "跳过RGI：$sample 已存在结果"
    fi
}

# 导出函数供parallel使用
export -f process_sample
export RGI_DB_PATH RGI_RESULT_DIR CONTIGS_DIR

# 并行处理所有样本
echo "开始并行处理 $SAMPLES ..."
parallel -j "$PARALLEL_JOBS" process_sample {} "$THREADS_PER_SAMPLE" ::: $SAMPLES

echo "🎉 所有样本RGI注释完成！"
echo "📁 RGI结果位置: $RGI_RESULT_DIR/*_arg.txt"