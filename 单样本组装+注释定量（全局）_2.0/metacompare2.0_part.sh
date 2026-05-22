#!/bin/bash
# 脚本名称：C_5_batch_metacompare2.sh
# 功能：批量使用 metacompare2.0 对所有单样本 contigs 进行抗性组风险评估

# ================= 用户配置区域 =================

# 设置 MetaCompare 2.0 脚本路径
METACOMPARE_SCRIPT="/database/work/zryan/software/MetaCompare2.0/metacompare.py"

# 设置单样本 contigs 输入目录
INPUT_DIR="./D_a1_SampleContigs/SContigs"

# 设置总体输出根目录
BASE_OUTPUT_DIR="./Risk_metacompare2"
mkdir -p "$BASE_OUTPUT_DIR"

# 设置线程数
THREADS=72 

# 设置风险评分模式：0 = 同时计算生态与人类健康风险
RISK_MODE=0

# ================= 新增修复逻辑 =================
# 定义数据库的绝对路径
DATABASE_DIR="/database/work/zryan/software/MetaCompare2.0/metacmpDB"

# 检查工作目录下是否已有软链接或文件夹，如果有则删除（防止报错）
if [ -e "./metacmpDB" ]; then
    rm -f "./metacmpDB"
fi

# 在当前工作目录创建软链接，指向真实的数据库目录
# 这一步是关键：MetaCompare 2.0 会去当前目录找 ./metacmpDB
ln -s "$DATABASE_DIR" "./metacmpDB"

echo "已建立数据库软链接: ./metacmpDB -> $DATABASE_DIR"
echo "----------------------------------------"

# ================= 脚本逻辑区域 =================
echo "开始批量 MetaCompare 2.0 分析..."
echo "输入目录: $INPUT_DIR"
echo "输出目录: $BASE_OUTPUT_DIR"

# 仅处理指定样本（最小修改：替换原 for 循环）
# 目标样本：T5_2,T5_3,T5_4,T5_5,T5_6,T6_1,T6_2,T6_3,T6_4,T6_5,T6_6
for SAMPLE_NAME in T5_2 T5_3 T5_4 T5_5 T5_6 T6_1 T6_2 T6_3 T6_4 T6_5 T6_6; do

    # 构建对应的 contigs 文件路径
    CONTIGS_FILE="$INPUT_DIR/${SAMPLE_NAME}.contigs.fa"
    
    # 检查文件是否存在（防止目录为空报错）
    if [ ! -f "$CONTIGS_FILE" ]; then
        echo "未找到文件: $CONTIGS_FILE"
        continue
    fi

    # 为当前样本创建专属输出目录
    SAMPLE_OUTPUT_DIR="$BASE_OUTPUT_DIR/$SAMPLE_NAME"
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    echo "----------------------------------------"
    echo "正在处理样本: $SAMPLE_NAME"
    echo "输入文件: $CONTIGS_FILE"
    echo "输出目录: $SAMPLE_OUTPUT_DIR"

    # 运行 MetaCompare 2.0
    # 注意：这里使用双引号包裹变量，防止路径中有空格导致报错
    python "$METACOMPARE_SCRIPT" \
      -c "$CONTIGS_FILE" \
      -t "$THREADS" \
      -b "$RISK_MODE" \
      -o "$SAMPLE_OUTPUT_DIR"

    # 检查上一条命令是否成功执行
    if [ $? -eq 0 ]; then
        echo "样本 $SAMPLE_NAME 分析完成！"
    else
        echo "样本 $SAMPLE_NAME 分析失败！请检查日志。"
    fi

done

echo "----------------------------------------"
echo "指定样本处理完毕！"
echo "结果保存在: $BASE_OUTPUT_DIR"