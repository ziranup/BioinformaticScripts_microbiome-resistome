#!/bin/bash

# 设置路径
INPUT_GENOMES="./D_b5_MAGs/filtered_mags"
OUTPUT_ROOT="./Tree_MAGs_gtdbtk"
ANNOTATION_DIR="${OUTPUT_ROOT}/annotation"
TREE_BAC_DIR="${OUTPUT_ROOT}/tree_bacteria"
TREE_ARC_DIR="${OUTPUT_ROOT}/tree_archaea"
USER_TREE_DIR="${OUTPUT_ROOT}/user_only_trees"
ANNO_TMP_DIR="${ANNOTATION_DIR}/tmp"
BAC_TREE_TMP_DIR="${TREE_BAC_DIR}/tmp"
ARC_TREE_TMP_DIR="${TREE_ARC_DIR}/tmp"

# 创建输出目录
mkdir -p "${ANNOTATION_DIR}" "${TREE_BAC_DIR}" "${TREE_ARC_DIR}" "${USER_TREE_DIR}" "${ANNO_TMP_DIR}" "${BAC_TREE_TMP_DIR}" "${ARC_TREE_TMP_DIR}"


# 判断是否已存在 classify 结果，若存在则跳过(断点调试)
CLASSIFY_DONE=false
if [ -s "${ANNOTATION_DIR}/bin.bac120.summary.tsv" ] || [ -s "${ANNOTATION_DIR}/bin.ar53.summary.tsv" ]; then
    echo "[INFO] Existing GTDB-Tk classification results found. Skipping classify_wf."
    CLASSIFY_DONE=true
fi

# Step 1: 物种注释 (classify_wf) —— 仅在未完成时运行
if [ "$CLASSIFY_DONE" = false ]; then
    echo "[INFO] Running GTDB-Tk classify_wf for taxonomic annotation..."
    gtdbtk classify_wf \
        --genome_dir "${INPUT_GENOMES}" \
        --out_dir "${ANNOTATION_DIR}" \
        --tmpdir "${ANNO_TMP_DIR}" \
        --extension fa \
        --prefix bin \
        --cpus 72 \
        --pplacer_cpus 72 \
        --scratch_dir "${ANNO_TMP_DIR}" \
        --keep_intermediates

    # 检查 classify 是否成功
    if [ ! -f "${ANNOTATION_DIR}/bin.bac120.summary.tsv" ] && [ ! -f "${ANNOTATION_DIR}/bin.ar53.summary.tsv" ]; then
        echo "[ERROR] classify_wf failed or no genomes detected."
        exit 1
    fi
else
    echo "[INFO] Proceeding with existing classification results."
fi

# Step 2: 判断是否存在细菌或古菌，分别建树
HAS_BAC=false
HAS_ARC=false

if [ -s "${ANNOTATION_DIR}/bin.bac120.summary.tsv" ]; then
    # 排除 header 行后是否有数据
    if [ $(wc -l < "${ANNOTATION_DIR}/bin.bac120.summary.tsv") -gt 1 ]; then
        HAS_BAC=true
    fi
fi

if [ -s "${ANNOTATION_DIR}/bin.ar53.summary.tsv" ]; then
    if [ $(wc -l < "${ANNOTATION_DIR}/bin.ar53.summary.tsv") -gt 1 ]; then
        HAS_ARC=true
    fi
fi

# Step 3: 细菌 de novo 树（如果需要全域，外群选古菌）
# 注意：可以在kraken中查找p__Patescibacteriota是否存在，如果不存在，就用这个作为外群（gtdbtk示例）
if [ "$HAS_BAC" = true ]; then
    echo "[INFO] Building bacterial de novo tree..."
    gtdbtk de_novo_wf \
        --genome_dir "${INPUT_GENOMES}" \
        --out_dir "${TREE_BAC_DIR}" \
        --tmpdir "${BAC_TREE_TMP_DIR}" \
        --extension fa \
        --bacteria \
        --outgroup_taxon p__Patescibacteriota \
        --prefix bin \
        --cpus 72 \
        --keep_intermediates

else
    echo "[INFO] No bacterial MAGs detected. Skipping bacterial tree."
fi

# Step 4: 古菌 de novo 树（如果需要全域，外群选细菌）
# 注意：可以在kraken中查找p__Altiarchaeota是否存在，如果不存在，就用这个作为外群（gtdbtk示例）
if [ "$HAS_ARC" = true ]; then
    echo "[INFO] Building archaeal de novo tree..."
    gtdbtk de_novo_wf \
        --genome_dir "${INPUT_GENOMES}" \
        --out_dir "${TREE_ARC_DIR}" \
        --tmpdir "${ARC_TREE_TMP_DIR}" \
        --extension fa \
        --archaea \
        --outgroup_taxon p__Altiarchaeota \
        --prefix bin \
        --cpus 72 \
        --keep_intermediates
else
    echo "[INFO] No archaeal MAGs detected. Skipping archaeal tree."
fi

# Step 5 (Optional): 仅用用户 MAGs 构建高精度树（基于 classify 产生的 MSA）
# Bacterial user-only tree
if [ "$HAS_BAC" = true ] && [ -f "${ANNOTATION_DIR}/align/bin.bac120.user_msa.fasta.gz" ]; then
    echo "[INFO] Inferring user-only bacterial tree from MSA..."
    gtdbtk infer \
        --msa_file "${ANNOTATION_DIR}/align/bin.bac120.user_msa.fasta.gz" \
        --out_dir "${USER_TREE_DIR}" \
        --cpus 72 \
        --prefix bin_bac_user
fi

# Archaeal user-only tree
if [ "$HAS_ARC" = true ] && [ -f "${ANNOTATION_DIR}/align/bin.ar53.user_msa.fasta.gz" ]; then
    echo "[INFO] Inferring user-only archaeal tree from MSA..."
    gtdbtk infer \
        --msa_file "${ANNOTATION_DIR}/align/bin.ar53.user_msa.fasta.gz" \
        --out_dir "${USER_TREE_DIR}" \
        --cpus 72 \
        --prefix bin_arc_user
fi

echo "[DONE] GTDB-Tk annotation and tree building completed."