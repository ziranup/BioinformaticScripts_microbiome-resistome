#!/bin/bash

# 脚本名称: run_arg_abundance.sh
# 功能描述: 在RGI完成后，单独计算ARGs丰度（CPM, TPM, FPKM）

# 设置进程数（每个样本分配12线程，并行6个样本，总72线程）
THREADS_PER_SAMPLE=12
PARALLEL_JOBS=6

# 设置路径
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"   # clean reads用于回帖（非压缩）
QC_READS_DIR="./D_1_QCreads"                 # QC_reads主要用于正确计算cpm（.gz 压缩）

# 输出目录
RGI_RESULT_DIR="./ARGs_CARD/rgi/result"
ABUNDANCE_RESULT_DIR="./ARGs_CARD/abundance/result"
ABUNDANCE_TMP_DIR="./ARGs_CARD/abundance/tmp"

# 创建输出目录
mkdir -p "$ABUNDANCE_RESULT_DIR" "$ABUNDANCE_TMP_DIR"

# === 检查 coverm 是否可用 ===
if ! command -v coverm &> /dev/null; then
    echo "错误: 'coverm' 命令未找到。"
    echo "请确保已激活包含 coverm 的 conda 环境。"
    exit 1
fi

# 步骤1: 清理临时目录并复制RGI结果
echo "清理临时目录并复制RGI结果..."
rm -rf "$ABUNDANCE_TMP_DIR"
mkdir -p "$ABUNDANCE_TMP_DIR"
cp "$RGI_RESULT_DIR"/*_arg.txt "$ABUNDANCE_TMP_DIR"/

# 提取所有样本名（从RGI结果文件推断）
# 示例：CK_1_arg.txt → CK_1
SAMPLES=$(ls "$ABUNDANCE_TMP_DIR"/*_arg.txt | xargs -n1 basename | sed 's/_arg\.txt$//')

# 定义处理单个样本的函数
process_sample() {
    local sample=$1
    local threads=$2

    # 路径定义
    local rgi_txt="$ABUNDANCE_TMP_DIR/${sample}_arg.txt"
    local r1_read="$CLEAN_READS_DIR/${sample}.hrm.1.fastq"
    local r2_read="$CLEAN_READS_DIR/${sample}.hrm.2.fastq"

    # CoverM临时与结果路径
    local coverm_count="$ABUNDANCE_TMP_DIR/${sample}_count.tsv"
    local coverm_tpm="$ABUNDANCE_TMP_DIR/${sample}_tpm.tsv"
    local coverm_fpkm="$ABUNDANCE_TMP_DIR/${sample}_fpkm.tsv"

    # 步骤2: 提取ARG基因ID（仅Strict/Perfect）- 使用正确的列索引
    local arg_gene_list="$ABUNDANCE_TMP_DIR/${sample}_arg_genes.txt"
    awk -F'\t' 'NR>1 && ($7 == "Strict" || $7 == "Perfect") {print $1}' "$rgi_txt" > "$arg_gene_list"

    # 调试：检查提取结果
    gene_count=$(wc -l < "$arg_gene_list" 2>/dev/null)
    echo "样本 $sample: 提取到 $gene_count 个ARG基因 (Cut_Off列检查)"

    # 若无ARG，则创建空丰度表
    if [ ! -s "$arg_gene_list" ] || [ "$gene_count" -eq 0 ]; then
        echo -e "ARG_Name\t${sample}_CPM\t${sample}_TPM\t${sample}_FPKM" > "$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"
        echo -e "No_ARG_Detected\t0.0\t0.0\t0.0" >> "$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"
        return 0
    fi

    # 步骤3: 从RGI结果提取Predicted_DNA构建ARG-specific reference
    # 关键修正：Predicted_DNA是第18列，Best_Hit_ARO是第9列
    local arg_orf_fasta="$ABUNDANCE_TMP_DIR/${sample}_arg_orfs.fasta"
    awk -F'\t' 'NR>1 && ($7 == "Strict" || $7 == "Perfect") {
        gsub(/ /, "_", $9);  
        print ">" $1 "|" $9 "\n" $18
    }' "$rgi_txt" | fold -w 60 > "$arg_orf_fasta"

    # 调试：检查FASTA文件
    if [ ! -s "$arg_orf_fasta" ]; then
        echo "警告: 样本 $sample 的ARG FASTA文件为空"
        echo -e "ARG_Name\t${sample}_CPM\t${sample}_TPM\t${sample}_FPKM" > "$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"
        echo -e "No_ARG_Detected\t0.0\t0.0\t0.0" >> "$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"
        return 0
    fi

    # 运行CoverM on ARG ORFs
    coverm contig \
        --reference "$arg_orf_fasta" \
        --coupled "$r1_read" "$r2_read" \
        --methods mean covered_fraction covered_bases variance count length reads_per_base rpkm tpm \
        --min-read-aligned-length 50 \
        --min-read-percent-identity 90 \
        --min-read-aligned-percent 80 \
        --threads "$threads" \
        --output-count "$coverm_count" \
        --output-tpm "$coverm_tpm" \
        --output-fpkm "$coverm_fpkm"

    # 步骤4: 合并为最终表格（每样本一个TSV）
    local final_table="$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"

    # 计算该样本的总 clean reads 数（双端）
    r1_clean="$QC_READS_DIR/${sample}.R1.qualified.fastq.gz"
    r2_clean="$QC_READS_DIR/${sample}.R2.qualified.fastq.gz"

    if [[ -f "$r1_clean" && -f "$r2_clean" ]]; then
        total_reads_r1=$(zcat "$r1_clean" | wc -l)
        total_reads_r2=$(zcat "$r2_clean" | wc -l)
        total_clean_reads=$(( (total_reads_r1 + total_reads_r2) / 4 ))
    else
        echo "⚠️ 警告: 样本 $sample 的 clean reads 缺失，尝试 fallback 到 D_2_HRreads..."
        r1_fallback="$CLEAN_READS_DIR/${sample}.hrm.1.fastq"
        r2_fallback="$CLEAN_READS_DIR/${sample}.hrm.2.fastq"
        if [[ -f "$r1_fallback" && -f "$r2_fallback" ]]; then
            total_clean_reads=$(( ($(wc -l < "$r1_fallback") + $(wc -l < "$r2_fallback")) / 4 ))
        else
            echo "❌ 错误: 无法获取样本 $sample 的 clean reads 总数！"
            total_clean_reads=1
        fi
    fi

    # 合并三列并计算CPM
    paste "$coverm_count" "$coverm_tpm" "$coverm_fpkm" | \
    awk -v total_clean="$total_clean_reads" -v sample="$sample" '
        NR==1 {
            print "ARG_Name\t" sample "_CPM\t" sample "_TPM\t" sample "_FPKM"
            next
        }
        {
            orf_id = $1
            split(orf_id, a, "\\|")
            arg_name = (length(a[2]) > 0) ? a[2] : a[1]
            count = $2; tpm = $4; fpkm = $6
            if (total_clean > 0) {
                cpm = (count / total_clean) * 1000000
            } else {
                cpm = 0.0
            }
            printf "%s\t%.6f\t%.6f\t%.6f\n", arg_name, cpm, tpm, fpkm
        }' > "$final_table"

    echo "✅ 样本 $sample 处理完成 (总 clean reads: $total_clean_reads)"
}

# 导出函数供parallel使用
export -f process_sample
export CLEAN_READS_DIR QC_READS_DIR ABUNDANCE_TMP_DIR ABUNDANCE_RESULT_DIR

# 并行处理所有样本
echo "开始并行处理 $SAMPLES ..."
parallel -j "$PARALLEL_JOBS" process_sample {} "$THREADS_PER_SAMPLE" ::: $SAMPLES

# 合并所有样本的丰度表
echo "正在合并所有样本的ARG丰度表..."
HEADER="ARG_Name"
for sample in $SAMPLES; do
    HEADER="${HEADER}\t${sample}_CPM\t${sample}_TPM\t${sample}_FPKM"
done
echo -e "$HEADER" > "$ABUNDANCE_RESULT_DIR/ALL_SAMPLES_ARG_abundance_combined.tsv"

# 获取所有ARG名称（去重）
ALL_ARG_NAMES=$(awk 'NR>1 {print $1}' "$ABUNDANCE_RESULT_DIR"/*_abundance.tsv | sort -u)

# 逐行构建合并表
while IFS= read -r arg; do
    line="$arg"
    for sample in $SAMPLES; do
        file="$ABUNDANCE_RESULT_DIR/${sample}_abundance.tsv"
        vals=$(awk -v a="$arg" '$1==a {print $2"\t"$3"\t"$4}' "$file")
        if [ -z "$vals" ]; then
            vals="0.0\t0.0\t0.0"
        fi
        line="$line\t$vals"
    done
    echo -e "$line"
done <<< "$ALL_ARG_NAMES" >> "$ABUNDANCE_RESULT_DIR/ALL_SAMPLES_ARG_abundance_combined.tsv"

echo "🎉 所有样本ARG丰度分析完成！"
echo "📁 单样本结果: $ABUNDANCE_RESULT_DIR/*_abundance.tsv"
echo "📊 合并表格: $ABUNDANCE_RESULT_DIR/ALL_SAMPLES_ARG_abundance_combined.tsv"