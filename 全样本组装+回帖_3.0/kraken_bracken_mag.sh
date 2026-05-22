#!/bin/bash

# 设置工作目录和数据库路径
MAGS_DIR="./D_b5_MAGs/filtered_mags"
OUTPUT_BASE="./Species_MAGs_kraken_bracken"
KRAKEN_DB="/database/work/zryan/biodb/kraken/Standard_260226"
KRAKEN_TOOLS_DIR="./Species_MAGs_kraken_bracken/KrakenTools"

# 创建输出子目录
mkdir -p "${OUTPUT_BASE}/kraken2"
mkdir -p "${OUTPUT_BASE}/bracken"
mkdir -p "${OUTPUT_BASE}/mpa"

# 定义单个MAG处理函数
process_mag() {
    local mag_file="$1"
    local base_name=$(basename "$mag_file" .fa)

    # Kraken2 注释
    kraken2 \
        --db "${KRAKEN_DB}" \
        --threads 72 \
        --report "${OUTPUT_BASE}/kraken2/${base_name}.kraken2.report" \
        --output "${OUTPUT_BASE}/kraken2/${base_name}.kraken2.result" \
        "$mag_file"

    # Bracken 校正（物种水平，read length=150）
    bracken \
        -d "${KRAKEN_DB}" \
        -i "${OUTPUT_BASE}/kraken2/${base_name}.kraken2.report" \
        -o "${OUTPUT_BASE}/bracken/${base_name}.kraken2.bracken" \
        -w "${OUTPUT_BASE}/bracken/${base_name}.kraken2.bracken.report" \
        -r 150 \
        -l S

    # 转换为 mpa 格式
    python3 "${KRAKEN_TOOLS_DIR}/kreport2mpa.py" \
        -r "${OUTPUT_BASE}/bracken/${base_name}.kraken2.bracken.report" \
        -o "${OUTPUT_BASE}/mpa/${base_name}.kraken2.bracken.mpa.report"
}

# 导出函数供 parallel 使用
export -f process_mag
export KRAKEN_DB
export OUTPUT_BASE
export KRAKEN_TOOLS_DIR


# 获取所有 MAG 文件列表
MAG_LIST=("${MAGS_DIR}"/*.fa)

# 使用 GNU parallel 并行处理（每个任务分配8线程，总核心数不超过72）
# 若系统无 parallel，可替换为普通 for 循环
printf '%s\n' "${MAG_LIST[@]}" | parallel -j 1 process_mag {}

# 合并所有 mpa 文件
python3 "${KRAKEN_TOOLS_DIR}/combine_mpa.py" \
    -i "${OUTPUT_BASE}/mpa/"*.mpa.report \
    -o "${OUTPUT_BASE}/combined_mpa.report"

echo "✅ 所有 MAG 物种注释完成！合并结果位于: ${OUTPUT_BASE}/combined_mpa.report"


















### 单个MAG模板

# # kraken2注释
# kraken2 --db /data/database/kraken2/gtdb/release202 --threads 50 /data/meta/metawrap/bin/XHMHH_bin.154.fa --report /data/meta/metawrap/bin/bin.154_kraken2gtdb.report --output /data/meta/metawrap/bin/bin.154_kraken2gtdb.result

# # Braken校正
# bracken -d /data/database/kraken2/gtdb/release202 -i /data/meta/metawrap/bin/bin.154_kraken2gtdb.report -o /data/meta/metawrap/bin/bin.154_kraken2gtdb.report.bracken -w /data/meta/metawrap/bin/bin.154_kraken2gtdb.report.bracken.report -r 150 -l S

# # Braken的report格式转换成--use-mpa-style格式
# kreport2mpa.py -r /data/meta/metawrap/bin/bin.154_kraken2gtdb.report.bracken.report -o /data/meta/metawrap/bin/bin.154_kraken2gtdb.report.bracken.mpa.report

# # 多个样本mpa文件合并
# combine_mpa.py -i /data/meta/metawrap/bin/*.mpa.report -o  /data/meta/metawrap/bin/combined_mpa.report

# # mpa文件格式转换（可选）
# transfer.pl /data/meta/metawrap/bin/combined_mpa.report > /data/meta/metawrap/bin/combined_mpa_transfer.txt