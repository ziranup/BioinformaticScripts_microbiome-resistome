#!/bin/bash
# 脚本名称：11_metabat2_binning_custom.sh
# 功能：基于统一 contigs 和各样本去宿主 reads 进行 MetaBAT2 分箱
# 输入：
#   - Contigs: ./D_3_contigs/mg_6.contigs.fa
#   - 去宿主 reads: ./D_2_HRreads/cleanreads/{sample}.hrm.1/2.fastq
# 输出：
#   - MAGs: ./D_4_bins/bin.*

set -e

# 路径定义
contigs_fa="./D_3_contigs/mg_6.contigs.fa"
reads_dir="./D_2_HRreads/cleanreads"
output_dir="./D_4_bins"
tmp_dir="./D_4_bins/tmp_dir"

# 创建输出与临时目录
mkdir -p "${output_dir}" "${tmp_dir}"

# 1. 构建 Bowtie2 索引（基于你的 contigs）
echo "🔍 正在构建 Bowtie2 索引..."
bowtie2-build "${contigs_fa}" "${tmp_dir}/contigs"

# 2. 定义单样本比对函数
process_sample() {
  local sample_id=$1
  local reads_dir=$2
  local tmp_dir=$3

  # 输入文件
  local r1="${reads_dir}/${sample_id}.hrm.1.fastq"
  local r2="${reads_dir}/${sample_id}.hrm.2.fastq"
  local sam_out="${tmp_dir}/${sample_id}.sam"
  local bam_out="${tmp_dir}/${sample_id}.bam"

  # 比对
  bowtie2 -x "${tmp_dir}/contigs" -1 "${r1}" -2 "${r2}" \
          -p 12 --very-sensitive -S "${sam_out}"

  # 转 BAM + 排序
  samtools view -@ 6 -bS "${sam_out}" | samtools sort -@ 6 -o "${bam_out}"

  # 删除 SAM 以节省空间
  rm -f "${sam_out}"
}

# 导出函数供 parallel 使用
export -f process_sample

# 3. 自动提取样本 ID 列表（从 .hrm.1.fastq 文件）
samples=$(ls "${reads_dir}"/*.hrm.1.fastq | xargs -n1 basename | sed 's/\.hrm\.1\.fastq$//')

# 4. 并行比对所有样本（使用 6 个并行任务，每个用 12 线程，总线程 ≤72）
echo "🧬 正在并行比对所有样本到 contigs..."
parallel -j 6 process_sample {} "${reads_dir}" "${tmp_dir}" ::: $samples

# 5. 合并 BAM 覆盖度生成 depth.txt
echo "📊 正在汇总覆盖度..."
jgi_summarize_bam_contig_depths --outputDepth "${tmp_dir}/depth.txt" "${tmp_dir}"/*.bam

# 6. 运行 MetaBAT2 分箱
echo "📦 正在运行 MetaBAT2 分箱..."
metabat2 -t 72 -i "${contigs_fa}" -a "${tmp_dir}/depth.txt" -o "${output_dir}/bin"

# 7. 清理临时 BAM 文件（可选，若磁盘紧张）
# rm -f "${tmp_dir}"/*.bam

echo "✅ MetaBAT2 分箱完成！"
echo "📁 MAGs 位于: ${output_dir}/"
echo "   文件名示例: bin.0.fa, bin.1.fa, ..."