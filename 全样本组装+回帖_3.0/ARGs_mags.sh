#!/bin/bash

# 数据库路径（已建索引）
SARG_DB="/database/work/zryan/biodb/SARG/SARG_3.2.4/SARG_protein.dmnd"

# 输入：ORF 蛋白质文件目录
PROTEINS_DIR="./D_b6_ORFs_MAGs"

# 输出：比对结果目录
OUTPUT_DIR="./ARGs_MAGs_SARG"
mkdir -p "$OUTPUT_DIR"

# 临时目录（确保有写权限）
TMPDIR="./tmp_diamond"
mkdir -p "$TMPDIR"

# 并行运行 DIAMOND
export SARG_DB OUTPUT_DIR TMPDIR

parallel -j 36 '
    protein_file={}
    sample_name=$(basename "$protein_file" _proteins.faa)
    
    diamond blastp \
        --db "$SARG_DB" \
        --query "$protein_file" \
        --out "$OUTPUT_DIR/${sample_name}_SARG.tsv" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
        --threads 2 \
        --evalue 1e-7 \
        --max-target-seqs 1 \
        --sensitive \
        --tmpdir "$TMPDIR" \
        --quiet \
        --header simple
' ::: "$PROTEINS_DIR"/*_proteins.faa

echo "DIAMOND 比对完成，结果位于 $OUTPUT_DIR"