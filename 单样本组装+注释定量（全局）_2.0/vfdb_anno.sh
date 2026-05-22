#!/bin/bash


VFDB_DB="/database/work/zryan/biodb/VFDB/VFDB_setB_pro.dmnd"
GENE_DIR="./D_a3_SampleNonredundantGenes"
OUT_DIR="./VFs_VFDB/VFDB_annotation"
THREADS=72


mkdir -p "$OUT_DIR"

annotate_sample() {
  local sample=$1
  local gene_file="$GENE_DIR/${sample}_nonredundant_genes.fna"
  local out_file="$OUT_DIR/${sample}.vfdb.tsv"

  diamond blastx \
    --db "$VFDB_DB" \
    --query "$gene_file" \
    --out "$out_file" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-6 \
    --threads "$THREADS" \
    --max-target-seqs 1 \
    --ultra-sensitive
}

# 获取样本名（安全方式）
mapfile -t samples < <(ls "$GENE_DIR"/*_nonredundant_genes.fna 2>/dev/null | xargs -n1 basename | sed 's/_nonredundant_genes.fna$//')

# 改用 for 循环（不再需要 export -f 或 parallel）
for sample in "${samples[@]}"; do
  echo "Processing $sample..."
  annotate_sample "$sample"
done

echo "VFDB 注释完成！结果在 $OUT_DIR"