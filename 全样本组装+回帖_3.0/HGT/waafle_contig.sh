#!/bin/bash

# 这个脚本用于contig级别的HGT注释
# 脚本有问题，contig级别请直接用命令行注释。

set -e

# 设置路径
INPUT_GENOMES="./D_b2_CoORFs/genes.fna"
OUTPUT_DIR="./HGTs_contig_waafle"
WAFFLE_DB_FILE="/database/work/zryan/biodb/WAAFLE/chocophlan.v202210_202403.waafledb/chocophlan.v202210_202403.waafledb"
WAFFLE_DB_TSV="/database/work/zryan/biodb/WAAFLE/chocophlan.v202210_202403.taxonomy.tsv"
MERGED_READS_1="./D_b1_CoContigs/mergereads/all.R1.fastq"
MERGED_READS_2="./D_b1_CoContigs/mergereads/all.R2.fastq"

mkdir -p $(dirname "${OUTPUT_FILE}")


# 第一步：Homology-based search with waafle_search
waafle_search "${INPUT_GENOMES}" "${WAFFLE_DB_FILE}" --out "${OUTPUT_FILE}/contigs.blastout" --threads 80 

# 第二步：Gene calling with waafle_genecaller
waafle_genecaller "${OUTPUT_DIR}/contigs.blastout"

# 第三步：Identify candidate LGT events with waafle_orgscorer
waafle_orgscorer \
  "${INPUT_GENOMES}" \
  "${OUTPUT_DIR}/contigs.blastout" \
  "${OUTPUT_DIR}/contigs.gff" \
  "${WAFFLE_DB_TSV}"

# 第四步：Filter out misassembled contigs
# waafle_junctions \
#   "${INPUT_GENOMES}" \
#   "${OUTPUT_DIR}/contigs.gff" \
#   --reads1 "${MERGED_READS_1}" \
#   --reads2 "${MERGED_READS_2}"
