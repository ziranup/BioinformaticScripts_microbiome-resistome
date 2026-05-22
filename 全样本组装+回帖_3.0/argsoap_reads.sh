#!/bin/bash
# 脚本名称：4_argsoap_reads.sh
# 功能：基于去宿主 reads 预测抗性基因（ARGs）
# 前置数据准备：
#   - 去宿主 reads: sample1.hrm.1.fastq ... sample4.hrm.2.fastq
#   - SARG 数据库路径已配置（args_oap 内部调用）
# 输出：
#   - output_dir/stage1/{sample}
#   - output_dir/stage2/{sample}

# args_oap 并行抗性基因预测脚本（20251125）
# 输入：./2_bowtie/seq_hrm/xxx.hrm.1.fastq 与 xxx.hrm.2.fastq（共30样本）
# 输出：./4_argsoap/stage1_result/ 和 ./4_argsoap/stage2_result/
# 注意：args_oap 要求输入文件夹中 reads 命名为 sample_1.fastq / sample_2.fastq

set -e  # 遇错退出

# 1. 创建输出目录
mkdir -p ./4_argsoap/stage1_result
mkdir -p ./4_argsoap/stage2_result
mkdir -p ./4_argsoap/tmp_input  # 用于存放重命名后的 reads

# 2. 定义单样本处理函数
process_sample() {
    local sample=$1
    local input_dir="./2_bowtie/seq_hrm"
    local tmp_input="./4_argsoap/tmp_input/${sample}"
    local stage1_out="./4_argsoap/stage1_result/${sample}"
    local stage2_out="./4_argsoap/stage2_result/${sample}"
    local sarg_db="/database/work/zryan/biodb/SARG/Short_subdatabase"

    # 创建临时输入目录
    mkdir -p "$tmp_input"

    # 创建符号链接或复制（推荐硬链接或cp，避免路径含.的问题；这里用cp更稳妥）
    cp "${input_dir}/${sample}.hrm.1.fastq" "${tmp_input}/${sample}_1.fastq"
    cp "${input_dir}/${sample}.hrm.2.fastq" "${tmp_input}/${sample}_2.fastq"

    # Stage One: 比对到 SARG 数据库
    args_oap stage_one \
        -i "$tmp_input" \
        -o "$stage1_out" \
        -f fastq \
        -t 12

    # Stage Two: 功能注释与定量
    args_oap stage_two \
        -i "$stage1_out" \
        -o "$stage2_out" \
        -t 12

    # 可选：清理临时文件（节省空间）
    rm -rf "$tmp_input"
}

# 3. 导出函数供 parallel 使用
export -f process_sample

# 4. 提取样本名（从 .hrm.1.fastq 文件）
# 示例：ABC.hrm.1.fastq → ABC
samples=$(ls ./2_bowtie/seq_hrm/*.hrm.1.fastq | xargs -n1 basename | sed 's/\.hrm\.1\.fastq$//')

# 5. 并行处理（-j 6 表示同时跑6个样本，每个用12线程 → 总72线程）
echo "检测到以下样本："
echo "$samples"
echo "开始并行运行 args_oap (stage_one + stage_two) ..."

parallel -j 6 process_sample ::: $samples

echo "✅ 所有样本 ARGs 分析完成！"
echo "Stage1 结果：./4_argsoap/stage1_result/"
echo "Stage2 结果：./4_argsoap/stage2_result/"