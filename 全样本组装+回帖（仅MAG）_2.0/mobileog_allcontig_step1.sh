#!/bin/bash
# 脚本名称：MGEs_step1.sh
# 功能：使用 mobileOG-db 对ORF蛋白质序列进行 MGE 注释，ai建议用非冗余，但是微信推送上教程都用的是基于contig的all.faa(蛋白质序列)
# 注意：mobileOG-db 是蛋白质数据库，必须使用蛋白质序列进行比对

set -e

# 创建输出目录
mkdir -p ./MGEs_mobileOG/mobileOG/result
mkdir -p ./MGEs_mobileOG/mobileOG/tmp

# 检查输入文件类型
echo "检查输入文件..."
head -n 2 ./D_7_nonredundantGenes/nonredundant_genes.fasta

# 使用原始蛋白质文件进行注释
ORF_PROTEINS="./D_6_ORFs/proteins.faa"
MOBILEOG_DB="./MGEs_mobileOG/mobileOG-db_1.6/mobileOG-db_beatrix-1.6.All.dmnd"
OUTPUT_FILE="./MGEs_mobileOG/mobileOG/result/mobileog_hits.tsv"

echo "开始使用 DIAMOND blastp 对蛋白质序列进行 MGE 注释..."
echo "使用的输入文件: $ORF_PROTEINS"

# 验证输入文件确实是蛋白质序列（更可靠的检测方法）
# 提取前10行的序列部分（跳过以>开头的标题行），检查是否包含典型的氨基酸字母
if head -n 20 "$ORF_PROTEINS" | grep -v '^>' | tr -d '\n' | grep -q '[MRNDHK]' ; then
    echo "✓ 检测到蛋白质序列（包含典型氨基酸字母 M, R, N, D, H, K 等）"
else
    echo "⚠ 警告：输入文件可能不是蛋白质序列！"
    echo "显示前几行序列内容以供检查："
    head -n 10 "$ORF_PROTEINS" | grep -v '^>'
    echo ""
    echo "请确认您使用的是蛋白质序列文件（.faa 后缀）"
    exit 1
fi

# 运行 DIAMOND blastp 比对
diamond blastp \
    --db "$MOBILEOG_DB" \
    --query "$ORF_PROTEINS" \
    --out "$OUTPUT_FILE" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
    --threads 72 \
    --evalue 1e-6 \
    --max-target-seqs 1 \
    --header simple \
    --tmpdir ./MGEs_mobileOG/mobileOG/tmp \
    --sensitive \
    --verbose

echo "✅ MGE 注释完成！结果保存至 $OUTPUT_FILE"
echo "接下来请运行丰度计算脚本：MGEs_step2.py"