#!/bin/bash
# 单样本MEGAHIT组装脚本（20260502）
# 目标：对30个样本分别组装contigs，使用meta-sensitive模式

# 1. 创建主输出目录和临时目录
mkdir -p ./D_8_SampleContigs
mkdir -p ./D_8_SampleContigs/tmp

# 2. 定义单样本组装函数
assemble_sample() {
  local sample=$1
  local r1="./D_2_HRreads/cleanreads/${sample}.hrm.1.fastq"
  local r2="./D_2_HRreads/cleanreads/${sample}.hrm.2.fastq"
  local out_dir="./D_8_SampleContigs/${sample}"
  local tmp_dir="./D_8_SampleContigs/tmp/${sample}"

  # 检查输入文件是否存在
  if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
    echo "⚠️ 跳过样本 $sample：R1或R2文件缺失"
    return 1
  fi

  # 👇 只创建 tmp_dir（MEGAHIT 不会自动创建 --tmp-dir）
  mkdir -p "$tmp_dir"
  # ❌ 不要创建 out_dir！让 MEGAHIT 自己创建

  echo "🧬 开始组装样本: $sample"

  megahit \
    -1 "$r1" \
    -2 "$r2" \
    --presets meta-sensitive \
    -m 0.7 \
    --mem-flag 1 \
    -t 12 \
    --out-dir "$out_dir" \
    --out-prefix "$sample" \
    --keep-tmp-files \
    --tmp-dir "$tmp_dir"

  # 验证输出
  if [ -f "${out_dir}/${sample}.contigs.fa" ]; then
    echo "✅ 样本 $sample 组装成功！"
  else
    echo "❌ 样本 $sample 组装失败！"
  fi
}



# 3. 导出函数供parallel使用
export -f assemble_sample

# 4. 提取样本名列表（去重）
samples=$(ls ./D_2_HRreads/cleanreads/*.hrm.1.fastq | xargs -n1 basename | sed 's/\.hrm\.1\.fastq$//' | sort -u)

# 5. 并行运行（-j 6 表示同时运行6个样本，每个占12线程，共72核）
parallel -j 6 assemble_sample ::: $samples

echo "🎉 所有单样本MEGAHIT组装完成！结果位于 ./D_8_SampleContigs/"