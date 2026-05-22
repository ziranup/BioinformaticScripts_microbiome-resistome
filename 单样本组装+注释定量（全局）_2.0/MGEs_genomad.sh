#!/bin/bash
# genomad_MGEs_prediction.sh
# Author: YZR
# Date: 2026-05-04

# 1. 创建输出与临时目录
mkdir -p ./MGEs_genomad
mkdir -p ./MGEs_genomad/tmp

# 2. 定义处理单个样本的函数
run_genomad() {
  local sample_id=$1
  local input_fasta="./D_a1_SampleContigs/SContigs/${sample_id}.contigs.fa"
  local output_dir="./MGEs_genomad/${sample_id}"
  local db_path="/database/work/zryan/biodb/geNomad/genomad_db"


  # 运行 genomad end-to-end
  genomad end-to-end \
    --restart \
    --quiet \
    --splits 4 \
    --threads 80 \
    --enable-score-calibration \
    --max-fdr 0.05 \
    "$input_fasta" \
    "$output_dir" \
    "$db_path"
}


# 3. 导出函数供 parallel 使用
export -f run_genomad

# 4. 获取样本ID列表（去除 .contigs.fa 后缀）
samples=$(ls ./D_a1_SampleContigs/SContigs/*.contigs.fa | xargs -n1 basename | sed 's/\.contigs\.fa$//')

# 5. 并行执行
parallel -j 1 run_genomad ::: $samples

echo "✅ 所有样本的 MGEs 预测已完成！结果位于 ./MGEs_genomad/"