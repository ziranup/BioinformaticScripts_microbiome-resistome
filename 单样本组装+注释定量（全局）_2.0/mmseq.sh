#!/bin/bash

# 获取 mmseq.sh 脚本自身所在的绝对目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 基于脚本目录，构建数据库的绝对路径
DB="${SCRIPT_DIR}/metacmpDB/GTDB/gtdb"

# 执行 mmseqs 命令
mmseqs createdb "$1" "${2}sample.contigs"
mmseqs taxonomy "${2}sample.contigs" "$DB" "${2}sample.assignments" "${2}sample.tmpFolder" --tax-lineage 1 --majority 0.7 --vote-mode 1 --lca-mode 3 --orf-filter 1 --split-memory-limit 50G
mmseqs createtsv "${2}sample.contigs" "${2}sample.assignments" "${2}${3}"

# 清理临时文件 (取消注释以启用)
# rm "${2}sample.contigs"*
# rm "${2}sample.assignments"*
# rm -r "${2}sample.tmpFolder"