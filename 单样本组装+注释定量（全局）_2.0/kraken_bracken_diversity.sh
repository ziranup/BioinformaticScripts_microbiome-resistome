#!/bin/bash
# Kraken2 + Bracken + Alpha Diversity pipeline (20260430)
# Input: ./D_2_HRreads/cleanreads/*.hrm.[12].fastq
# Output:
#   - Kraken reports: ./Community_kraken_bracken/kraken/
#   - Bracken outputs: ./Community_kraken_bracken/bracken/
#   - Alpha diversity: ./Community_kraken_bracken/diversity/

# ============= 1. 创建输出目录 =============
mkdir -p ./Community_kraken_bracken/kraken
mkdir -p ./Community_kraken_bracken/bracken
mkdir -p ./Community_kraken_bracken/diversity

# ============= 2. 定义 Kraken2 函数 =============
process_kraken() {
  local sample=$1
  local r1_in="./D_2_HRreads/cleanreads/${sample}.hrm.1.fastq"
  local r2_in="./D_2_HRreads/cleanreads/${sample}.hrm.2.fastq"
  local report_out="./Community_kraken_bracken/kraken/${sample}.kreport"

  kraken2 --threads 72 \
          --db /database/work/zryan/biodb/kraken/Standard_260226 \
          --paired "$r1_in" "$r2_in" \
          --report "$report_out"
        #   --memory-mapping
}
export -f process_kraken

# ============= 3. 提取样本名 =============
samples=$(ls ./D_2_HRreads/cleanreads/*.hrm.1.fastq | xargs -n1 basename | sed 's/\.hrm\.1\.fastq//')

# ============= 4. 并行运行 Kraken2 =============
echo "Running Kraken2 on all samples..."
parallel -j 1 process_kraken ::: $samples




# ============= 5. 定义 Bracken 函数 =============
process_bracken() {
  local sample=$1
  local kraken_in="./Community_kraken_bracken/kraken/${sample}.kreport"
  local bracken_out="./Community_kraken_bracken/bracken/${sample}.S.bracken"
  local bracken_report="./Community_kraken_bracken/bracken/${sample}.S.bracken.kreport"

  bracken -d /database/work/zryan/biodb/kraken/gtdb_v226_250609 \
          -i "$kraken_in" \
          -o "$bracken_out" \
          -w "$bracken_report" \
          -l S \
          -t 72
}
export -f process_bracken

# ============= 6. 并行运行 Bracken =============
echo "Running Bracken on all Kraken reports..."
parallel -j 1 process_bracken ::: $samples

# ============= 7. 定义 Alpha 多样性函数 =============
# 注意：alpha_diversity.py 每次只能计算一个指标，因此需循环调用8次
process_alpha_all() {
  local sample=$1
  local bracken_in="./Community_kraken_bracken/bracken/${sample}.S.bracken"
  local out_dir="./Community_kraken_bracken/diversity"

  # List of alpha indices to compute
  indices=("chao1" "ace" "shannon" "simpson" "richness" "pielou_e" "invsimpson" "observed_species")

  for idx in "${indices[@]}"; do
    python ./Community_kraken_bracken/KrakenTools/DiversityTools/alpha_diversity.py \
           -f "$bracken_in" \
           -a "$idx" > "${out_dir}/${sample}.alpha_${idx}.txt"
  done
}
export -f process_alpha_all

# ============= 8. 并行计算所有 Alpha 多样性 =============
echo "Computing 8 alpha diversity indices for all samples..."
parallel -j 30 process_alpha_all ::: $samples

# ============= 9. 合并所有 Alpha 多样性为一个表格 =============
# 每个样本生成8个文件，现在合并成一个宽表：行=样本，列=8个指标
python3 <<'EOF'
import os
import pandas as pd

div_dir = "./Community_kraken_bracken/diversity"
indices = ["chao1", "ace", "shannon", "simpson", "richness", "pielou_e", "invsimpson", "observed_species"]

samples = set()
for f in os.listdir(div_dir):
    if f.endswith(".txt") and any(idx in f for idx in indices):
        sample = f.split(".alpha_")[0]
        samples.add(sample)

df_list = []
for sample in sorted(samples):
    row = {"Sample": sample}
    for idx in indices:
        file_path = os.path.join(div_dir, f"{sample}.alpha_{idx}.txt")
        if os.path.exists(file_path):
            with open(file_path) as fh:
                val = fh.read().strip()
                try:
                    row[idx.capitalize()] = float(val)
                except:
                    row[idx.capitalize()] = val
        else:
            row[idx.capitalize()] = None
    df_list.append(row)

df = pd.DataFrame(df_list)
df.to_csv("./Community_kraken_bracken/alpha_diversity_combined.csv", index=False)
print("Combined alpha diversity table saved to ./Community_kraken_bracken/alpha_diversity_combined.csv")
EOF

echo "✅ Pipeline completed! Final alpha diversity table: ./Community_kraken_bracken/alpha_diversity_combined.csv"