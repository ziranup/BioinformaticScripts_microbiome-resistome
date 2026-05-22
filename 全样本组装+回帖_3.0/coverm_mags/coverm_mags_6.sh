#!/bin/bash

# 设置路径
MAGS_DIR="./D_b8_Coverage_MAGs/MAGs_6"
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
OUTPUT_DIR="./D_b8_Coverage_MAGs"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 明确列出所有样本名（与你提供的文件一致）
SAMPLES=(
  CK_1 CK_2 CK_3 CK_4 CK_5 CK_6
  T1_1 T1_2 T1_3 T1_4 T1_5 T1_6
  T3_1 T3_2 T3_3 T3_4 T3_5 T3_6
  T5_1 T5_2 T5_3 T5_4 T5_5 T5_6
  T6_1 T6_2 T6_3 T6_4 T6_5 T6_6
)

# 构建 reads 参数列表（空格分隔的 R1 R2 R1 R2 ...）
COVERM_READS=""
for sample in "${SAMPLES[@]}"; do
  R1="${CLEAN_READS_DIR}/${sample}.hrm.1.fastq"
  R2="${CLEAN_READS_DIR}/${sample}.hrm.2.fastq"
  # 检查文件是否存在（可选，增强鲁棒性）
  if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
    echo "Warning: Missing reads for sample ${sample}"
  fi
  COVERM_READS+="${R1} ${R2} "
done

# 获取所有 MAG 文件，按文件号排序（确保 bin.10 在 bin.2 之后）
for mag_file in $(ls "${MAGS_DIR}"/bin.*.fa | sort -V); do
  # 提取基础名，如 bin.100.fa → bin.100
  mag_name=$(basename "${mag_file}" .fa)
  output_file="${OUTPUT_DIR}/${mag_name}_mag_coverage.tsv"
  
  echo "Processing MAG: ${mag_name}"

  # 设置临时目录（避免 /tmp 空间不足）
  tmp_dir="${OUTPUT_DIR}/tmp_${mag_name}"
  mkdir -p "${tmp_dir}"

  # 运行 CoverM
  TMPDIR="${tmp_dir}" \
  coverm contig \
    --reference "${mag_file}" \
    --coupled \
    ${COVERM_READS} \
    --methods mean covered_fraction covered_bases variance count length reads_per_base rpkm tpm \
    --min-read-aligned-length 50 \
    --min-read-percent-identity 90 \
    --min-read-aligned-percent 80 \
    --threads 80 \
    --output-file "${output_file}"

  # 清理临时目录
  rm -rf "${tmp_dir}"

  echo "Done: ${output_file}"
done

echo "All MAGs processed. Results saved in ${OUTPUT_DIR}"