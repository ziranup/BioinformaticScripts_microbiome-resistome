#!/bin/bash

# 脚本名称: run_arg_analysis.sh
# 功能描述: 基于RGI和CARD数据库进行ARGs注释，并计算其在多个样本中的丰度（CPM, TPM, RPKM/FPKM）。

# 注意：contig模式和protein模式不是以编码识别的。蛋白质的genes.fna也应当以protein模式识别，因为这回跳过识别ORFs的步骤。




# 设置进程数
THREADS=72

# 设置输入文件路径
ALL_SAMPLE_CONTIGS="./D_3_contigs/mg_6.contigs.fa"
ALL_SAMPLE_ORFS_DIR="./D_6_ORFs"
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
NONREDUNDANT_GENE_FILE="/data/yanziran/mg_6_new/D_7_nonredundantGenes/nonredundant_genes.fasta"

# 设置数据库路径
RGI_DB_PATH="/database/work/zryan/biodb/CARD/card.json"

# 设置输出目录
RGI_RESULT_DIR="./ARGs_CARD/rgi/result"
RGI_TMP_DIR="./ARGs_CARD/rgi/tmp"
ABUNDANCE_RESULT_DIR="./ARGs_CARD/abundance/result"
ABUNDANCE_TMP_DIR="./ARGs_CARD/abundance/tmp"

# 创建所有必要的输出目录
mkdir -p "$RGI_RESULT_DIR" "$RGI_TMP_DIR" "$ABUNDANCE_RESULT_DIR" "$ABUNDANCE_TMP_DIR"

echo "==================== 步骤1: RGI注释ARGs ===================="

# 导入rgi数据库
rgi load --card_json $RGI_DB_PATH --local

# 使用RGI对非冗余基因集进行注释。
rgi main \
    -i "$NONREDUNDANT_GENE_FILE" \
    -o "$RGI_RESULT_DIR/nonredundant_genes_arg" \
    -t contig \
    -a DIAMOND \
    -n $THREADS \
    --include_nudge \
    --local \
    --debug

echo "RGI注释完成，结果保存在: $RGI_RESULT_DIR"

echo "==================== 步骤2: 使用CoverM计算基因丰度 ===================="

# === 新增：检查 coverm 是否可用 ===
if ! command -v coverm &> /dev/null; then
    echo "错误: 'coverm' 命令未找到。"
    echo "请确保您已激活包含 coverm 的 conda 环境。"
    echo "安装命令: mamba install -c conda-forge -c bioconda coverm"
    exit 1
fi
# === 检查结束 ===

# 创建 CoverM 所需的临时子目录
mkdir -p "${ABUNDANCE_TMP_DIR}/allmethods"
mkdir -p "${ABUNDANCE_TMP_DIR}/metabat"
mkdir -p "${ABUNDANCE_TMP_DIR}/coverage_histogram"

# 准备动态的 reads 列表
READS_LIST=""
for read1_file in "$CLEAN_READS_DIR"/*hrm.1.fastq; do
    sample_name=$(basename "$read1_file" | sed 's/\.hrm\.1\.fastq$//')
    read2_file="$CLEAN_READS_DIR/${sample_name}.hrm.2.fastq"
    READS_LIST="${READS_LIST} ${read1_file} ${read2_file}"
done

# 运行coverm,其他所有的方法
TMPDIR="${ABUNDANCE_TMP_DIR}/allmethods" coverm contig \
  --reference "${NONREDUNDANT_GENE_FILE}" \
  --coupled \
  ${READS_LIST} \
  --methods mean covered_fraction covered_bases variance count length reads_per_base rpkm tpm \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${ABUNDANCE_RESULT_DIR}/STDOUT" 


# # 运行coverm,coverage_histogram方法，因为不能与其他一起运行
# TMPDIR="${ABUNDANCE_TMP_DIR}/coverage_histogram" coverm contig \
#   --reference "${NONREDUNDANT_GENE_FILE}" \
#   --coupled \
#   ${READS_LIST} \
#   --methods coverage_histogram \
#   --min-read-aligned-length 50 \
#   --min-read-percent-identity 90 \
#   --min-read-aligned-percent 80 \
#   --threads 72 \
#   --output-file "${ABUNDANCE_RESULT_DIR}/STDOUT_coverage_histogram" 


# # 运行coverm,单独运行metabat方法.因为不能与其他一起运行
# TMPDIR="${ABUNDANCE_TMP_DIR}/metabat" coverm contig \
#   --reference "${NONREDUNDANT_GENE_FILE}" \
#   --coupled \
#   ${READS_LIST} \
#   --methods metabat \
#   --min-read-aligned-length 50 \
#   --min-read-percent-identity 90 \
#   --min-read-aligned-percent 80 \
#   --threads 72 \
#   --output-file "${ABUNDANCE_RESULT_DIR}/STDOUT_metabat" 


echo "✅ CoverM 基因丰度定量完成！结果在 ${ABUNDANCE_RESULT_DIR}/"
echo "💡 输出包含以下指标：mean, coverage_histogram, covered_fraction, covered_bases, variance, length, count, reads_per_base, rpkm, tpm, metabat"

# 说明：参数设置原因
# --min-read-aligned-length 50 \
# 宏基因组 reads 通常为 150 bp（PE），50 bp 是可靠比对的下限。过低会引入噪声；过高会丢失短 contig 的覆盖信息。类似研究常用 50–75 bp（Sczyrba et al., 2017, *Nature Methods*）。
# --min-read-percent-identity 90 \
# 宏基因组中物种多样性高，90–95% identity 是区分近缘菌株的常用阈值。90% 平衡了灵敏度与特异性，避免将 reads 错配到远缘同源序列（Nayfach et al., 2016, *Genome Research*）。
# --min-read-aligned-percent 80 \
# 	要求至少 80% 的 read 被比对上，防止部分比对（partial alignment）导致的假阳性覆盖。这在短 contigs（<1 kb）中尤为重要（Bishara et al., 2018, *Nature Biotechnology*）。

echo "==================== 步骤3: 汇总ARG丰度 (CPM, TPM, RPKM) ===================="
# 1. 从RGI结果中提取ARG基因ID列表
echo "正在从RGI结果中提取ARG信息..."
ARG_INFO_FILE="$RGI_TMP_DIR/arg_info.txt"
awk -F'\t' 'NR>1 && ($7 == "Strict" || $7 == "Perfect") {print $1"\t"$8}' \
    "$RGI_RESULT_DIR/nonredundant_genes_arg.txt" > "$ARG_INFO_FILE"

# 2. 读取CoverM的count文件来计算CPM，并整合TPM和RPKM
COVERM_COUNT_FILE="$ABUNDANCE_RESULT_DIR/STDOUT_count.tsv"
COVERM_TPM_FILE="$ABUNDANCE_RESULT_DIR/STDOUT_tpm.tsv"
COVERM_RPKM_FILE="$ABUNDANCE_RESULT_DIR/STDOUT_rpkm.tsv"
FINAL_ARG_TABLE="$ABUNDANCE_RESULT_DIR/ARG_abundance_table_with_units.tsv"

# 首先，创建一个包含所有ARG的初始表格（基于count）
TEMP_ARG_COUNT="$ABUNDANCE_TMP_DIR/arg_count.tsv"
awk -v arg_file="$ARG_INFO_FILE" '
BEGIN {
    while ((getline line < arg_file) > 0) {
        split(line, a, "\t");
        arg_map[a[1]] = a[2];
    }
    close(arg_file);
}
NR==1 { 
    print "ARG_Name\t" $0 > "'$TEMP_ARG_COUNT'";
    next;
}
$1 in arg_map {
    print arg_map[$1] "\t" $0 >> "'$TEMP_ARG_COUNT'";
}' "$COVERM_COUNT_FILE"

# 计算每个样本的总计数（用于CPM）
TOTAL_COUNTS="$ABUNDANCE_TMP_DIR/total_counts.txt"
tail -n +2 "$COVERM_COUNT_FILE" | awk '{for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) printf "%s\n", sum[i]}' > "$TOTAL_COUNTS"

# 生成最终表格，包含 CPM, TPM, RPKM
# 表头
{
    printf "ARG_Name"
    for (i=2; i<=NF; i++) {
        sample = $i;
        printf "\t%s_CPM\t%s_TPM\t%s_RPKM", sample, sample, sample;
    }
    printf "\n";
} > "$FINAL_ARG_TABLE"

# 处理数据行
tail -n +2 "$TEMP_ARG_COUNT" | while IFS=$'\t' read -r arg_name counts; do
    # 将counts字符串转换为数组
    IFS=' ' read -ra COUNT_ARRAY <<< "$counts"
    
    printf "%s" "$arg_name"
    for i in "${!COUNT_ARRAY[@]}"; do
        if [ $i -eq 0 ]; then continue; fi # 跳过第一个元素（基因名）
        count_val=${COUNT_ARRAY[$i]}
        sample_col=$((i+1)) # 在total_counts中对应的列
        
        # 获取总计数
        total_count=$(sed -n "${sample_col}p" "$TOTAL_COUNTS")
        
        # 计算CPM
        if [ "$total_count" -gt 0 ]; then
            cpm=$(awk "BEGIN {printf \"%.6f\", ($count_val / $total_count) * 1000000}")
        else
            cpm="0.000000"
        fi
        
        # 从TPM和RPKM文件中提取对应值
        tpm_val=$(awk -v gene="$arg_name" -v col="$sample_col" 'NR>1 && $1==gene {print $col}' "$COVERM_TPM_FILE")
        rpkm_val=$(awk -v gene="$arg_name" -v col="$sample_col" 'NR>1 && $1==gene {print $col}' "$COVERM_RPKM_FILE")
        
        printf "\t%s\t%s\t%s" "$cpm" "$tpm_val" "$rpkm_val"
    done
    printf "\n"
done >> "$FINAL_ARG_TABLE"

echo "ARG丰度汇总完成，最终结果保存在: $FINAL_ARG_TABLE"
echo "==================== 分析全部完成 ===================="