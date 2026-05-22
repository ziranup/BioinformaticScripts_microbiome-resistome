#!/bin/bash
# quantify_vf_with_coverm.sh

HR_READS_DIR="./D_2_HRreads/cleanreads"
GENE_DIR="./D_a3_SampleNonredundantGenes"
ANNOT_DIR="./VFs_VFDB/VFDB_annotation"
COVERM_OUT="./VFs_VFDB/VFDB_abundance"
THREADS_PER_SAMPLE=20

mkdir -p "$COVERM_OUT"

quantify_sample() {
  local sample=$1
  local r1="$HR_READS_DIR/${sample}.hrm.1.fastq"
  local r2="$HR_READS_DIR/${sample}.hrm.2.fastq"
  local ref="$GENE_DIR/${sample}_nonredundant_genes.fna"
  local annot="$ANNOT_DIR/${sample}.vfdb.tsv"
  local out_prefix="$COVERM_OUT/${sample}"

  # Step 1: 用 CoverM 计算每个基因的 coverage & read count
  coverm contig \
    --methods mean covered_bases variance count length reads_per_base rpkm tpm \
    --input-fasta "$ref" \
    --coupled "$r1" "$r2" \
    --min-read-aligned-length 50 \
    --min-read-percent-identity 90 \
    --min-read-aligned-percent 80 \
    --output-format sparse \
    --threads "$THREADS_PER_SAMPLE" \
    --output-file "$out_prefix.coverm.tsv"

  # Step 2: 提取 VF 基因（根据注释结果中的 qseqid）
  awk '{print $1}' "$annot" | sort | uniq > "$out_prefix.vf_ids.txt"
  join -t $'\t' <(sort "$out_prefix.vf_ids.txt") <(sort "$out_prefix.coverm.tsv") > "$out_prefix.vf_abundance.tsv"
}

export -f quantify_sample

samples=$(ls "$ANNOT_DIR"/*.vfdb.tsv | xargs -n1 basename | sed 's/.vfdb.tsv//')

parallel -j 4 quantify_sample ::: $samples

echo "CoverM 定量完成！结果在 $COVERM_OUT"