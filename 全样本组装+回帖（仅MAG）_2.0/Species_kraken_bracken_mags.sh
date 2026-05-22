#!/bin/bash
# 使用 Kraken2 + Bracken 注释 MAGs 并计算其在原始样本中的丰度
# 输入:
#   - MAGs: ./D_b5_MAGs/filtered_mags/*.fa
#   - Clean reads: ./D_2_HRreads/cleanreads/*.hrm.[12].fastq
# 数据库:
#   - GTDB Kraken2 DB: /database/work/zryan/biodb/kraken/gtdb_v226_250609
# 输出:
#   - ./MAGs_annotation/kraken/        : MAGs 的 Kraken2 报告
#   - ./MAGs_annotation/bracken/       : Bracken 物种级重估计结果
#   - ./MAGs_annotation/species_info.tsv : 每个 MAG 的主导物种信息
#   - ./MAGs_annotation/abundance/     : 每个样本中各 MAG 物种的 Bracken 丰度表

# ============= 1. 显式定义并导出路径 =============
MAGS_DIR="./D_b5_MAGs/filtered_mags"
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
KRAKEN_DB="/database/work/zryan/biodb/kraken/gtdb_v226_250609"
OUTPUT_DIR="./MAGs_kraken_bracken"

# 👇 关键：导出为环境变量，供 parallel 子进程使用
export KRAKEN_DB

KRANKEN_OUT="$OUTPUT_DIR/kraken"
BRACKEN_OUT="$OUTPUT_DIR/bracken"
ABUNDANCE_OUT="$OUTPUT_DIR/abundance"
mkdir -p "$KRANKEN_OUT" "$BRACKEN_OUT" "$ABUNDANCE_OUT"

# ============= 2. 定义函数：注释单个 MAG =============
annotate_single_mag() {
  local mag_fasta="$1"
  local basename=$(basename "$mag_fasta")
  local mag_id="${basename%.fa}"

  local kraken_report="$KRANKEN_OUT/${mag_id}.kreport"
  local bracken_out="$BRACKEN_OUT/${mag_id}.S.bracken"
  local bracken_report="$BRACKEN_OUT/${mag_id}.S.bracken.kreport"

  # 使用 Kraken2 对 MAG (fasta) 进行分类
  kraken2 --threads 72 \
          --db "$KRAKEN_DB" \
          --report "$kraken_report" \
          "$mag_fasta"

  # 使用 Bracken，指定读长为 150（匹配你的数据）
  bracken -d "$KRAKEN_DB" \
          -i "$kraken_report" \
          -o "$bracken_out" \
          -w "$bracken_report" \
          -l S \
          -t 72
}
export -f annotate_single_mag

# ============= 3. 获取所有 MAG 文件 =============
mag_files=("$MAGS_DIR"/*.fa)
if [ ! -e "${mag_files[0]}" ]; then
  echo "❌ 错误：未在 $MAGS_DIR 中找到任何 .fa 文件！"
  exit 1
fi

# ============= 4. 串行注释所有 MAGs（替换 parallel 以避免内存溢出） =============
echo "🔍 共检测到 ${#mag_files[@]} 个 MAG，开始串行注释（避免内存溢出）..."
for mag in "${mag_files[@]}"; do
  echo "Processing: $(basename "$mag")"
  annotate_single_mag "$mag"
done

# ============= 5. 提取每个 MAG 的主导物种信息 =============
echo "🧾 正在提取 MAG 的主导物种信息..."
python3 << EOF
import os
import pandas as pd

bracken_dir = "$BRACKEN_OUT"
output_file = "$OUTPUT_DIR/species_info.tsv"

records = []
for f in os.listdir(bracken_dir):
    if f.endswith(".S.bracken.kreport"):
        mag_id = f.replace(".S.bracken.kreport", "")
        report_path = os.path.join(bracken_dir, f)
        try:
            df = pd.read_csv(report_path, sep='\t', header=None, comment='#')
            df.columns = ['percent', 'reads', 'tax_reads', 'kmers', 'tax_kmers', 'taxid', 'rank', 'name']
            # 优先取 species，其次 genus
            species_row = df[df['rank'] == 'species']
            if not species_row.empty:
                top = species_row.iloc[0]
                assigned = top['name'].strip()
                rank = 'species'
            else:
                genus_row = df[df['rank'] == 'genus']
                if not genus_row.empty:
                    top = genus_row.iloc[0]
                    assigned = f"{top['name'].strip()} (genus)"
                    rank = 'genus'
                else:
                    assigned = "Unassigned"
                    rank = 'unassigned'
            records.append({
                'MAG': mag_id,
                'Assigned_Taxon': assigned,
                'Rank': rank,
                'Confidence_Percent': float(top['percent']) if 'top' in locals() else 0.0
            })
        except Exception as e:
            print(f"⚠️  跳过 {f}: {e}")

pd.DataFrame(records).to_csv(output_file, sep='\t', index=False)
print(f"✅ 物种注释表已保存至: {output_file}")
EOF

# ============= 6. 为每个样本运行 Kraken2 + Bracken（复用你已有的流程，但仅生成报告）============
# 注意：这里我们不重新运行整个 pipeline，而是假设你已有 Community_kraken_bracken/bracken/ 下的 .S.bracken 文件
# 如果尚未运行，请先运行你提供的 reads 注释流程！

# 检查 reads 注释结果是否存在
if [ ! -d "./Community_kraken_bracken/bracken" ]; then
  echo "⚠️  警告：未检测到 ./Community_kraken_bracken/bracken/ 目录。"
  echo "   请先运行你提供的 reads 注释流程，再继续此脚本！"
  exit 1
fi

# ============= 7. 构建“MAG 物种 -> 样本丰度”矩阵 =============
# 思路：每个 MAG 被赋予一个物种名（如 "Pseudomonas sp001234567"），我们在每个样本的 Bracken 物种表中查找该名称的丰度

echo "📊 正在构建 MAG 物种在各样本中的丰度矩阵..."

python3 << EOF
import os
import pandas as pd
import glob

# 加载 MAG 物种分配表
species_df = pd.read_csv("$OUTPUT_DIR/species_info.tsv", sep='\t')
# 只保留成功注释到 species 或 genus 的 MAG
species_df = species_df[species_df['Rank'].isin(['species', 'genus'])].copy()

# 提取纯物种名（去除 "(genus)" 后缀用于匹配）
species_df['Search_Name'] = species_df['Assigned_Taxon'].str.replace(r' $$(genus)$$', '', regex=True)

# 获取所有样本的 Bracken 物种文件
bracken_files = glob.glob("./Community_kraken_bracken/bracken/*.S.bracken")
if not bracken_files:
    raise FileNotFoundError("未找到任何 .S.bracken 文件！请确保 reads 注释已完成。")

# 构建样本名列表
samples = []
for f in bracken_files:
    sample_name = os.path.basename(f).replace(".S.bracken", "")
    samples.append(sample_name)

# 初始化丰度矩阵
abundance_matrix = pd.DataFrame(index=species_df['MAG'], columns=samples, data=0.0)

# 逐个读取样本的 Bracken 结果
for f in bracken_files:
    sample_name = os.path.basename(f).replace(".S.bracken", "")
    df = pd.read_csv(f, sep='\t')
    # Bracken 列: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
    taxon_to_frac = dict(zip(df['name'], df['fraction_total_reads']))
    
    # 为每个 MAG 查找其物种在该样本中的丰度
    for idx, row in species_df.iterrows():
        taxon = row['Search_Name']
        if taxon in taxon_to_frac:
            abundance_matrix.loc[row['MAG'], sample_name] = taxon_to_frac[taxon]
        # 否则保持为 0.0

# 保存结果
abundance_matrix.to_csv("$ABUNDANCE_OUT/MAG_species_abundance_matrix.tsv", sep='\t')
print(f"✅ MAG 物种丰度矩阵已保存至: $ABUNDANCE_OUT/MAG_species_abundance_matrix.tsv")
print("   行 = MAG ID, 列 = 样本, 值 = 该 MAG 所代表物种在样本中的相对丰度（0~1）")
EOF

echo "✅ 全流程完成！"
echo "📁 主要输出目录: $OUTPUT_DIR"
echo "   - 物种注释: $OUTPUT_DIR/species_info.tsv"
echo "   - 丰度矩阵: $ABUNDANCE_OUT/MAG_species_abundance_matrix.tsv"