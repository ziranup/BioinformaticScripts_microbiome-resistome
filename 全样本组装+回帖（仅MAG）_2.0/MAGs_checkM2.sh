#!/bin/bash
# 脚本名称：19_checkm2_filter_mags.sh
# 功能：使用 CheckM2 评估 MAG 质量，并筛选高质量 MAGs（Completeness ≥70%, Contamination ≤10%）
# 输入：
#   - 原始 MAGs 目录：/database/work/zryan/analysis/mg_6_new/D_5_MAGs/
#   - MAG 文件命名格式：bin.*.fa （如 bin.1.fa, bin.23.fa）
# 输出（统一在 ./D_9_filtered_MAGs/ 下）：
#   - checkm_output/        → CheckM2 原始输出（含 quality_report.tsv）
#   - filtered_mags/        → 符合质量标准的 MAGs (.fa)
#   - checkm_filtered_quality.tsv → 筛选后的质量报告

set -e

# === 配置区 ===
INPUT_MAG_DIR="./D_4_bins"
OUTPUT_BASE="./D_5_MAGs"
# 数据库：需要定位到文件而不是所在路径
CHECKM2_DB_PATH="/database/work/zryan/biodb/checkm2_1.1.0/CheckM2_database/uniref100.KO.1.dmnd"

# 创建输出目录
mkdir -p "${OUTPUT_BASE}/checkm_output" "${OUTPUT_BASE}/filtered_mags"

# === Step 1: 运行 CheckM2 质量预测 ===
echo "🔬 正在运行 CheckM2 质量评估..."
checkm2 predict \
  --threads 72 \
  --input "${INPUT_MAG_DIR}" \
  --output-directory "${OUTPUT_BASE}/checkm_output" \
  --database_path "${CHECKM2_DB_PATH}" \
  --extension ".fa" \
  --force

# === Step 2: 筛选高质量 MAGs ===
QUALITY_REPORT="${OUTPUT_BASE}/checkm_output/quality_report.tsv"
FILTERED_TSV="${OUTPUT_BASE}/checkm_filtered_quality.tsv"

# 列顺序（根据 CheckM2 v1.1.0 输出）：
# Column 1: Name
# Column 2: Completeness
# Column 3: Contamination
# 所以条件：$2 >= 50 && $3 <= 10
awk -F'\t' 'NR==1 || ($2 >= 50 && $3 <= 10)' "${QUALITY_REPORT}" > "${FILTERED_TSV}"

# === Step 3: 复制合格 MAGs 到 filtered_mags/ ===
echo "📦 正在复制高质量 MAGs..."
while IFS= read -r bin_id; do
  if [ -n "$bin_id" ] && [ "$bin_id" != "Name" ]; then
    SRC_FILE="${INPUT_MAG_DIR}/${bin_id}.fa"
    if [ -f "$SRC_FILE" ]; then
      cp "$SRC_FILE" "${OUTPUT_BASE}/filtered_mags/"
    else
      echo "⚠️ 警告：未找到 ${SRC_FILE}，跳过。"
    fi
  fi
done < <(tail -n +2 "${FILTERED_TSV}" | cut -f1)

# === 完成提示 ===
TOTAL_MGS=$(wc -l < "${QUALITY_REPORT}" | awk '{print $1-1}')
FILTERED_COUNT=$(wc -l < "${FILTERED_TSV}" | awk '{print $1-1}')

echo ""
echo "✅ CheckM2 质控完成！"
echo "  - 总 MAGs 数量: ${TOTAL_MGS}"
echo "  - 高质量 MAGs (C≥70%, Cont≤10%): ${FILTERED_COUNT}"
echo "  - 质量报告: ${FILTERED_TSV}"
echo "  - 高质量 MAGs 已保存至: ${OUTPUT_BASE}/filtered_mags/"