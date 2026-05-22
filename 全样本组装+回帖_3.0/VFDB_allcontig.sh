#!/bin/bash
# 毒力因子（VF）注释与丰度计算流程（基于VFDB + CoverM）
# 作者：YZR
# 时间：2026-05-07

# =============== 配置路径 ===============
RAW_READS_DIR="./D_0_rawdata"
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
CONTIGS_FA="./D_b1_CoContigs/mg_6.contigs.fa"
PROTEINS_FA="./D_b2_CoORFs/proteins.faa"
VFDB_DMND="/database/work/zryan/biodb/VFDB/VFDB_setB_pro.dmnd"
OUTPUT_ROOT="./VFs_VFDB_allcontig"

# 创建输出目录
mkdir -p "${OUTPUT_ROOT}/1_vfdb_diamond"
mkdir -p "${OUTPUT_ROOT}/2_contig_to_VF"
mkdir -p "${OUTPUT_ROOT}/3_coverm_contig_abundance"
mkdir -p "${OUTPUT_ROOT}/4_VF_abundance_tables"

# =============== Step 1: DIAMOND 比对 ORF 到 VFDB ===============
echo "【Step 1】正在运行 DIAMOND 比对..."

diamond blastp \
  --db "${VFDB_DMND}" \
  --query "${PROTEINS_FA}" \
  --out "${OUTPUT_ROOT}/1_vfdb_diamond/proteins_vs_VFDB.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --evalue 1e-7 \
  --max-target-seqs 1 \
  --threads 72 \
  --sensitive

# =============== Step 2: 提取 contig → VF 映射 ===============
echo "【Step 2】正在提取 contig 到 VF 的映射关系..."

# Prodigal ORF 格式: contigID_xxx → 提取 contigID（去掉最后一个下划线及之后）
awk '{
    n = split($1, a, "_");
    if (n < 2) { print $1 "\t" $2; next }
    contig = a[1];
    for (i = 2; i <= n-1; i++) {
        contig = contig "_" a[i];
    }
    print contig "\t" $2
}' "${OUTPUT_ROOT}/1_vfdb_diamond/proteins_vs_VFDB.tsv" > "${OUTPUT_ROOT}/2_contig_to_VF/contig_to_VF.tsv"

# =============== Step 3: 使用 CoverM 计算 contig 覆盖度 ===============
echo "【Step 3】正在使用 CoverM 计算 contig 覆盖度..."

# 构建 CoverM reads 参数列表（按你提供的格式）
COVERM_READS=""
for sample in CK_1 CK_2 CK_3 CK_4 CK_5 CK_6 \
              T1_1 T1_2 T1_3 T1_4 T1_5 T1_6 \
              T3_1 T3_2 T3_3 T3_4 T3_5 T3_6 \
              T5_1 T5_2 T5_3 T5_4 T5_5 T5_6 \
              T6_1 T6_2 T6_3 T6_4 T6_5 T6_6; do
  R1="${CLEAN_READS_DIR}/${sample}.hrm.1.fastq"
  R2="${CLEAN_READS_DIR}/${sample}.hrm.2.fastq"
  COVERM_READS="${COVERM_READS} ${R1} ${R2}"
done

# 运行 CoverM（严格按你提供的格式）
TMPDIR="${OUTPUT_ROOT}/3_coverm_contig_abundance" coverm contig \
  --reference "${CONTIGS_FA}" \
  --coupled \
  ${COVERM_READS} \
  --methods mean covered_fraction covered_bases variance count length reads_per_base rpkm tpm \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_ROOT}/3_coverm_contig_abundance/coverm_contig_coverage.tsv" \
  --verbose

# =============== Step 4: 将 contig-level 丰度汇总到 VF-level ===============
echo "【Step 4】正在汇总 contig 丰度到 VF 丰度..."

# 准备临时文件
TMP_DIR=$(mktemp -d)
cp "${OUTPUT_ROOT}/3_coverm_contig_abundance/coverm_contig_coverage.tsv" "${TMP_DIR}/coverage.tsv"
cp "${OUTPUT_ROOT}/2_contig_to_VF/contig_to_VF.tsv" "${TMP_DIR}/map.tsv"

# 添加 header 到 map 文件
sed -i '1i\gene_id\tVF_id' "${TMP_DIR}/map.tsv"

# 修正后的 R 脚本
RSCRIPT="${OUTPUT_ROOT}/4_VF_abundance_tables/aggregate_vf_abundance.R"
cat > "${RSCRIPT}" << 'EOF'
# 加载必要库
library(reshape2)
library(dplyr)
library(tidyr)

# 读取数据
coverage <- read.delim("coverage.tsv", row.names=1, check.names=FALSE)
map <- read.delim("map.tsv", stringsAsFactors=FALSE)

# 转换为 data.frame 并添加 contig 列
cov_df <- data.frame(contig = rownames(coverage), coverage, row.names = NULL)

# 获取列名
cols <- colnames(cov_df)

# 分离 TPM 和 RPKM 列
tpm_cols <- cols[grepl("\\.tpm$", cols)]
rpkm_cols <- cols[grepl("\\.rpkm$", cols)]

# 处理 TPM 表
if(length(tpm_cols) > 0) {
  tpm_data <- cov_df[, c("contig", tpm_cols)]
  tpm_long <- melt(tpm_data, id.vars = "contig", variable.name = "sample_metric", value.name = "TPM")
  tpm_long$sample <- gsub("\\.tpm$", "", tpm_long$sample_metric)
  
  tpm_merged <- merge(tpm_long, map, by.x = "contig", by.y = "gene_id", all.x = FALSE)
  vf_tpm <- tpm_merged %>%
    group_by(sample, VF_id) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = sample, values_from = TPM, values_fill = 0)
  write.table(vf_tpm, "sample_VF_TPM.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}

# 处理 RPKM 表
if(length(rpkm_cols) > 0) {
  rpkm_data <- cov_df[, c("contig", rpkm_cols)]
  rpkm_long <- melt(rpkm_data, id.vars = "contig", variable.name = "sample_metric", value.name = "RPKM")
  rpkm_long$sample <- gsub("\\.rpkm$", "", rpkm_long$sample_metric)
  
  rpkm_merged <- merge(rpkm_long, map, by.x = "contig", by.y = "gene_id", all.x = FALSE)
  vf_rpkm <- rpkm_merged %>%
    group_by(sample, VF_id) %>%
    summarise(RPKM = sum(RPKM, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = sample, values_from = RPKM, values_fill = 0)
  write.table(vf_rpkm, "sample_VF_RPKM.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}
EOF

# 运行 R 脚本
cd "${TMP_DIR}"
Rscript "${RSCRIPT}"

# 移动结果
mv sample_VF_TPM.tsv "${OUTPUT_ROOT}/4_VF_abundance_tables/" 2>/dev/null || echo "警告：未生成 TPM 表（可能无匹配）"
mv sample_VF_RPKM.tsv "${OUTPUT_ROOT}/4_VF_abundance_tables/" 2>/dev/null || echo "警告：未生成 RPKM 表（可能无匹配）"

# 清理临时文件
rm -rf "${TMP_DIR}"

echo "✅ 所有步骤完成！"
echo "VF 丰度表位于：${OUTPUT_ROOT}/4_VF_abundance_tables/"