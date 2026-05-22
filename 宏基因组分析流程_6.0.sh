##########
# 宏基因组分析流程_5.0 #
##########
# 5.0版本增加了大量的流程。实际上不需要跑那么多，视情况而定。


############目录############
前言：服务器配置查看
第一部分：数据准备模块
## 子任务27：使用KneadData工具，输入raw reads进行质控+去除宿主操作
## 子任务1：下机序列判断（子任务27下位替代的步骤1）
## 子任务2：使用fastp对下机序列进行质控（子任务27下位替代的步骤2）
## 子任务3：去除宿主序列(bowtie2)
## 子任务4：宏基因组组装（使用megahit生成contigs）
## 子任务5：宏基因组分箱binning（使用MetaBAT2生成bins）
## 子任务6：使用prodigal得到ORFs（CDS）
## 子任务7：构建非冗余基因目录（用于丰度定量）
## 子任务8：将MAGs转化为MAG水平蛋白质序列
## 子任务9：过滤低质量 MAGs（checkM2）（可选，但是可能丢失信息）
## 子任务30：使用dRep工具，将MAGs进行去冗余，得到物种水平的代表性基因组（SRGs）
第二部分：基于reads层面的分析
## 子任务29：使用AGS-and-ACN-tools计算平均基因组大小（AGS）和16S rRNA基因平均拷贝数（ACN）
## 子任务10：使用kraken和bracken进行微生物分类和丰度定量
## 子任务11：ARGs_OAP 抗性基因注释（reads层面）
## 子任务12：HUMAnN3功能注释
第三部分：基于contigs、ORF层面的分析
## 子任务45：使用Quast工具，对基因组组装结果进行质量评估
## 子任务13：reads水平基因丰度定量（样本的回帖）
## 子任务14：毒力因子预测
## 子任务50：使用GO数据库，创建中间映射；并以另一个独立脚本给出GO中间映射到注释的操作流程
## 子任务15：KEGG通路注释
## 子任务47：使用MetaCyC数据库进行代谢通路注释
## 子任务31：使用mobileOG-db工具，将contigs级别的蛋白质序列（ORFs）进行MGEs注释
## 子任务46：使用TCDB数据库对微生物组转运功能基因进行注释
## 子任务40：使用CAZyme数据库进行碳降解基因注释
## 子任务39：使用NCycDB进行氮循环功能基因注释
## 子任务43：使用PCycDB数据库，对磷循环功能基因进行注释
## 子任务44：使用MCycDB数据库，分析微生物组的甲烷循环过程
## 子任务48：使用Swiss-Prot数据库，对蛋白质进行注释
## 子任务42：使用Pfam数据库对蛋白质进行结构域层面的注释
## 子任务49：使用NR数据库，对蛋白质进行注释
## 子任务41：使用InterProScan对蛋白进行注释和功能分析
## 子任务51：使用OrthoDB数据库对蛋白质进行注释
## 子任务16：自定义数据库的比对解释（归类到第三部分其他独立数据库注释子任务之后）
## 子任务52：使用metacompare2.0进行抗性组风险评估
第四部分：基于MAGs层面的分析
## 子任务17：MAGs丰度定量
## 子任务18/子任务34：使用GTDB-tk工具，对MAGs进行物种注释与建树
## 子任务19：可移动遗传原件MGEs功能注释（基于MAGs）
## 子任务20：使用 DIAMOND 映射 MAGs 蛋白质到 ARGs 数据库（如 SARG）
## 子任务24/子任务33：使用CarveMe进行微生物群落基因组尺度代谢建模（与33相同）
## 子任务25：群体基因组学（within-species variation）
## 子任务21：水平基因转移（HGT）检测
## 子任务35：使用DRAM工具进行微生物组的代谢基因注释（重新写脚本，使得其在一个脚本中，既进行ORFs层面的注释分析，又进行MAGs层面的注释分析）
第六部分：独立完整的生物信息学流程
## 子任务32:使用metaGEM工具，从宏基因组层面重建基因组尺度的代谢模型
## 子任务37：使用MetaWRAP工具，从宏基因组数据挖掘单菌基因组bins（从质控到分箱+注释的独立整体流程）
第七部分：病毒组专题分析
## 子任务22：病毒组挖掘（Virome）
## 子任务23：CRISPER-Cas系统识别
## 子任务36：使用DRAM-v工具进行病毒组代谢基因注释
第八部分：多组学联合分析
## 子任务26：宏基因组与宏转录组整合分析
第九部分：扩增子联合分析
## 子任务28：使用Barrnap工具，从宏基因组中提取出16sRNA
## 子任务38：SortMeRNA（V4）工具：宏基因组提取rRNA序列
############目录############









# 前言：服务器配置查看 ####

## 一、查看逻辑CPU数目（最大线程）####
grep -c ^processor /proc/cpuinfo
maxcore
# 我已经用alias把maxcore写进.bashrc了，很方便。

## 二、查看最大内存（最大运行参数）####
free -h
# 最容易导致报错的点，需要保守设置

## 三、查看分区剩余空间####
df -h

# 《分析思路：宏基因组典型流程》#### 
# Raw reads → (fastp&去宿主) → Clean reads
# Clean reads → (MEGAHIT / metaSPAdes) → Contigs（或 Scaffolds）
# Contigs/Scaffolds → (Prodigal) → genes.fasta（CDS）
# genes.fasta → (DIAMOND/eggNOG-mapper等) → 功能注释



第一部分：数据准备模块
## 子任务1：下机序列判断
# 这一步需要做的工作有：
# 验证md5是否与给的md5一致。
# 判断文件格式，是压缩的序列还是没有压缩的序列。
# 弄清楚分组
# 重命名测序的文件。通常不同的测序文件按照样本和双端序列分开。
Md5sum ./rawdata/*


## 子任务27：使用KneadData工具，输入raw reads进行质控+去除宿主操作
# 放在fastq-bowtie2之前，作为这两个 新增子任务的上位替代

#!/bin/bash
# 脚本名称：1_kneaddata_qc.sh
# 功能：使用 KneadData 对 raw reads 进行质控与去宿主
# 输入：
#   - ./raw_reads/sample1.R1.raw.fastq.gz, ./raw_reads/sample1.R2.raw.fastq.gz ...
# 宿主参考基因组索引路径：
#   - /database/host_db/hg38/bowtie2_index/hg38
# 输出目录：./1_kneaddata/
# 输出文件：sample1_humann.fastq, sample1_humann_1.fastq, sample1_humann_2.fastq ...
# 笔记：
# （1）目的：对原始宏基因组测序数据进行质量控制并去除宿主（如人类）DNA污染。这是一个完整的单步骤流程（包含Trimmomatic质控 + Bowtie2去宿主）。
# （2）适用性：完全适用于宏基因组学、微生物组宏基因组、肠道微生物宏基因组研究（尤其常见于人源样本）。
# （3）归类：第一部分：数据准备模块
# （4）工具与环境：主要工具为 kneaddata（依赖 Trimmomatic + Bowtie2）。可通过 conda 独立环境在 Linux 中完成。
# （5）数据准备：raw reads（sample*.R1.raw.fastq.gz, sample*.R2.raw.fastq.gz）；需提前构建宿主参考基因组索引（如 human hg38）。

set -e

# 创建输出目录
mkdir -p ./1_kneaddata

# 定义宿主数据库路径（需提前用 bowtie2-build 构建）
HOST_DB="/database/host_db/hg38/bowtie2_index/hg38"

# 遍历所有样本
for r1 in ./raw_reads/*R1.raw.fastq.gz; do
    r2=${r1/R1.raw.fastq.gz/R2.raw.fastq.gz}
    sample=$(basename "$r1" | sed 's/.R1.raw.fastq.gz//')
    
    echo "Processing $sample ..."
    
    kneaddata \
        --input "$r1" \
        --input "$r2" \
        --reference-db "$HOST_DB" \
        --output ./1_kneaddata \
        --threads 16 \
        --trimmomatic /opt/anaconda3/envs/kneaddata/share/trimmomatic-*/ \
        --bowtie2-options "--very-sensitive --dovetail" \
        --output-prefix "${sample}_hrm"
done

echo "✅ KneadData 完成！clean reads 保存至 ./1_kneaddata/"



## 子任务2：使用fastp对下机序列进行质控（子任务27下位替代的步骤1）

#!/bin/bash
# 脚本名称：1_fastp_qc.sh
# 功能：对原始双端 reads 进行质控
# 前置数据准备：
#   - 当前目录存在 sample1.R1.fastq, sample1.R2.fastq, ..., sample4.R1.fastq, sample4.R2.fastq
# 输出：
#   - output_dir/sample{1..4}.R1.qualified.fastq.gz
#   - output_dir/sample{1..4}.R2.qualified.fastq.gz
#   - output_dir/report/

# fastp双端处理脚本（20251112）（YZR)
# 笔记：
# parallel对.符号敏感，需要把fastp封装到函数中。GNU Parallel 对命令模板中的特殊字符（比如 .）解析敏感。当命令模板中包含 .（常见于文件路径或命名中，比如 ./rawdata、.R1.raw.fastq），Parallel 可能会误认为是特殊语法而报错。报错内容：parallel: Error: Command cannot contain the character . Use a function for that.
# 2.0笔记：如何节省空间
# fastp 支持读取.fastq.gz压缩文件（无需提前gunzip），也支持直接输出.fastq.gz压缩文件（通过--compression参数），且 bowtie2 确实可以直接读取.fastq.gz格式文件。因此只需对脚本做少量修改即可满足需求，核心改动为调整输入输出文件后缀和添加压缩参数，代码结构基本不变。
# 1. 创建输出目录（若不存在）
mkdir -p ./1_QC
mkdir -p ./1_QC/report
# 2. 定义处理单个样本的函数（核心：封装fastp命令）
process_sample() {
  local sample=$1  # 接收样本名参数
  # 输入文件路径（改动1：输入改为.fastq.gz）
  local r1_in="./rawdata/${sample}.R1.raw.fastq.gz"
  local r2_in="./rawdata/${sample}.R2.raw.fastq.gz"
  # 输出文件路径（改动2：输出改为.fastq.gz）
  local r1_out="./1_QC/${sample}.R1.qualified.fastq.gz"
  local r2_out="./1_QC/${sample}.R2.qualified.fastq.gz"
  # 报告路径
  local html_report="./1_QC/report/${sample}_fastp.html"
  local json_report="./1_QC/report/${sample}_fastp.json"
  # 运行fastp
  fastp -i "$r1_in" -I "$r2_in" \
        -o "$r1_out" -O "$r2_out" \
        -h "$html_report" -j "$json_report" \
        -w 12 \
        --compression 6 \  # 6为平衡压缩率和速度的常用值，可根据需求调整
        --verbose
}

# 3. 导出函数（让parallel能识别）
export -f process_sample
# 4. 提取样本名列表（改动4：适配.gz后缀的输入文件）
samples=$(ls ./rawdata/*.R1.raw.fastq.gz | xargs -n1 basename | sed 's/\.R1\.raw\.fastq\.gz//')
# 5. 并行调用函数处理所有样本（-j 并行数，根据服务器核心数调整）
parallel -j 6 process_sample ::: $samples
echo "所有样本质控完成！结果在./1_QC，报告在./1_QC/report"



## 子任务3：去除宿主序列(bowtie2)（子任务27下位替代的步骤2）
# 较新的方式，已经使用kneddata进行。kneddata仍然使用的是bowtie2的核心算法。但是是更加新的工具。

#!/bin/bash
# 脚本名称：2_bowtie2_remove_host.sh
# 功能：去除弹尾虫宿主序列（参考基因组：Folsomia candida）
# 前置数据准备：
#   - 已有质控后 reads: sample1.R1.qualified.fastq ... sample4.R2.qualified.fastq
#   - 当前目录存在宿主参考基因组 host_genome.fa（即 F. candida 基因组）
# 输出：
#   - output_dir/sample{1..4}.hrm.1.fastq
#   - output_dir/sample{1..4}.hrm.2.fastq

# bowtie2批量去除宿主序列脚本（20251115优化空间占用）
# 处理对象：./1_QC下的30对质控后样本（xxx.R1.qualified.fastq / xxx.R2.qualified.fastq）
# 输出路径：./2_bowtie/seq_hrm（去宿主序列）、./2_bowtie/sam_output（比对结果）
# 该脚本还需要用dos2unix进行处理。具体处理方法见Debugnote_dos2unix.sh
# 笔记：解决文件命名混乱的问题。
# 修改前输出文件名格式异常（如CK_1.hrm.R.1.fastq），是因为 --un-conc-gz 参数的前缀设置不当。
# 具体来说：bowtie2 的 --un-conc-gz 会在指定的前缀后自动添加 .1 和 .2 来区分双端文件（R1 和 R2）。
# 不如从这一步开始，不用R1和R2，直接用hrm.1.fastq.gz和hrm.2.fastq.gz
# 存在的问题：
# 目前的脚本输出格式仍然为.fastq，该脚本还需要用dos2unix进行处理。
# dos模式下的换行符仍然存在。
# 测试：
# 通过把生成的fastq文件下载下来，T3_2.hrm.1.fastq。在windows系统上，用记事本打开，如果不能够打开，则为实际为.gz的压缩包文件。如果能够打开，则说明下面的脚本参数无效。
# 结果：生成的T3_2.hrm.1.fastq文件实际上是一个.gz的文件，使用winrar打开，里面的T3_2.hrm.1文件实际上才是那个fastq文件。
# 解决（未测试）：
# 直接在--un-conc-gz参数后面补上.gz后缀。

# 1. 创建输出目录（若不存在）
mkdir -p ./2_bowtie/seq_hrm    # 存储去宿主后的fastq文件
mkdir -p ./2_bowtie/sam_output # 临时存储SAM文件（处理后自动删除）
# 2. 构建宿主基因组索引（仅首次运行时需要，若已构建可注释此步骤）
# 索引前缀为host_index，存储在./2_bowtie目录下
if [ ! -f "./2_bowtie/host_index.1.bt2" ]; then
  echo "开始构建宿主基因组索引..."
  bowtie2-build ./2_bowtie/GCF_fcandida.fna ./2_bowtie/host_index
  echo "宿主基因组索引构建完成！"
else
  echo "宿主基因组索引已存在，跳过构建步骤"
fi

# 3. 定义处理单个样本的函数（核心：封装bowtie2命令，压缩HRM+自动删除SAM）
process_sample() {
  local sample=$1  # 接收样本名（如CK_1）
  # 输入文件路径（质控后的序列）
  local r1_in="./1_QC/${sample}.R1.qualified.fastq"
  local r2_in="./1_QC/${sample}.R2.qualified.fastq"
# 输出文件路径（2.0）
  local sam_out="./2_bowtie/sam_output/${sample}_hrm.sam"  # 比对结果存放在临时SAM文件中
  # 新前缀：./2_bowtie/seq_hrm/CK_1.R（后续bowtie2会自动加1.hrm.fastq.gz和2.hrm.fastq.gz）
  local hrm_prefix="./2_bowtie/seq_hrm/${sample}"

# 运行bowtie2
  bowtie2 \
  --very-sensitive \
  -t \
  -p 12 \
  -x ./2_bowtie/host_index \
  -1 "$r1_in" \
  -2 "$r2_in" \
  -S "$sam_out" \
  --un-conc-gz "$hrm_prefix.hrm.fastq.gz"
  # 最终生成：${sample}.hrm.1.fastq.gz 和 ${sample}.hrm.2.fastq.gz

  # 处理完成后自动删除SAM文件（释放空间）
  rm -f "$sam_out"
  echo "样本${sample}处理完成：压缩去宿主序列在./2_bowtie/seq_hrm，SAM文件已自动删除"
}
# 4. 导出函数（让parallel能识别）
export -f process_sample
# 5. 提取样本名列表（从./1_QC下的R1文件中提取）
samples=$(ls ./1_QC/*.R1.qualified.fastq | xargs -n1 basename | sed 's/\.R1\.qualified\.fastq//')
# 6. 并行处理所有样本（-j 并行数，根据服务器总核心数调整，建议总线程数=服务器核心数）
# 例如：若服务器有48核，每个样本用24线程，则并行数=48/24=2，即-j 2
parallel -j 6 process_sample ::: $samples
echo "所有样本去宿主处理完成！"



## 子任务4：宏基因组组装（使用megahit生成contigs）

## 版本1：合并组装(不建议)
#!/bin/bash
# 03_megahit_assembly_all.sh
# 合并所有样本的 clean reads，仅用于后续 co-assembly
mkdir -p results/merged

time cat ./2_bowtie/seq_hrm/*.hrm.1.fastq > ./6_assembly/1_mergereads/all.R1.fastq
time cat ./2_bowtie/seq_hrm/*.hrm.2.fastq > ./6_assembly/1_mergereads/all.R2.fastq

echo "✅ Reads merged for co-assembly"

# 创建输出目录：不需要，不然会提示报错目录已存在。
# 创建临时文件输出目录：需要！
mkdir -p /database_new/work/zryan/TEMP_analysis/mg_6_megahit

# 检查输入文件是否存在
if [ ! -f ./6_assembly/1_mergereads/all.R1.fastq ] || [ ! -f ./6_assembly/1_mergereads/all.R2.fastq ]; then
  echo "❌ 输入文件缺失！"
  exit 1
fi

# 运行 MEGAHIT（专为大内存限制优化）
megahit \
  -1 ./6_assembly/1_mergereads/all.R1.fastq \
  -2 ./6_assembly/1_mergereads/all.R2.fastq \
  --presets meta-sensitive \
  -m 0.8 \
  --mem-flag 1 \
  -t 72 \
  --out-dir ./6_assembly/2_megahit \
  --out-prefix mg_6 \
  --keep-tmp-files \
  --tmp-dir /database_new/work/zryan/TEMP_analysis/mg_6_megahit

# 验证结果
if [ ! -f ./6_assembly/2_megahit/mg_6.contigs.fa ]; then
  echo "❌ MEGAHIT 组装失败！"
  exit 1
fi

echo "✅ MEGAHIT 完成！结果: ./6_assembly/2_megahit/mg_6.contigs.fa"


### 版本2：以样本输入组装（也不会保留样本信息）
#!/bin/bash
# 03_megahit_assembly_sample.sh

# 创建输出目录：不需要，不然会提示报错目录已存在。
# 创建临时文件输出目录：需要
mkdir -p /database_new/work/zryan/TEMP_analysis/mg_6_megahit

# 运行 MEGAHIT（高灵敏度,k-mer 范围建议从21开始）
megahit \
  -1 ./2_bowtie/seq_hrm/CK_1.hrm.1.fastq,./2_bowtie/seq_hrm/CK_2.hrm.1.fastq,./2_bowtie/seq_hrm/CK_3.hrm.1.fastq,./2_bowtie/seq_hrm/CK_4.hrm.1.fastq,./2_bowtie/seq_hrm/CK_5.hrm.1.fastq,./2_bowtie/seq_hrm/CK_6.hrm.1.fastq,./2_bowtie/seq_hrm/T1_1.hrm.1.fastq,./2_bowtie/seq_hrm/T1_2.hrm.1.fastq,./2_bowtie/seq_hrm/T1_3.hrm.1.fastq,./2_bowtie/seq_hrm/T1_4.hrm.1.fastq,./2_bowtie/seq_hrm/T1_5.hrm.1.fastq,./2_bowtie/seq_hrm/T1_6.hrm.1.fastq,./2_bowtie/seq_hrm/T3_1.hrm.1.fastq,./2_bowtie/seq_hrm/T3_2.hrm.1.fastq,./2_bowtie/seq_hrm/T3_3.hrm.1.fastq,./2_bowtie/seq_hrm/T3_4.hrm.1.fastq,./2_bowtie/seq_hrm/T3_5.hrm.1.fastq,./2_bowtie/seq_hrm/T3_6.hrm.1.fastq,./2_bowtie/seq_hrm/T5_1.hrm.1.fastq,./2_bowtie/seq_hrm/T5_2.hrm.1.fastq,./2_bowtie/seq_hrm/T5_3.hrm.1.fastq,./2_bowtie/seq_hrm/T5_4.hrm.1.fastq,./2_bowtie/seq_hrm/T5_5.hrm.1.fastq,./2_bowtie/seq_hrm/T5_6.hrm.1.fastq,./2_bowtie/seq_hrm/T6_1.hrm.1.fastq,./2_bowtie/seq_hrm/T6_2.hrm.1.fastq,./2_bowtie/seq_hrm/T6_3.hrm.1.fastq,./2_bowtie/seq_hrm/T6_4.hrm.1.fastq,./2_bowtie/seq_hrm/T6_5.hrm.1.fastq,./2_bowtie/seq_hrm/T6_6.hrm.1.fastq \
  -2 ./2_bowtie/seq_hrm/CK_1.hrm.2.fastq,./2_bowtie/seq_hrm/CK_2.hrm.2.fastq,./2_bowtie/seq_hrm/CK_3.hrm.2.fastq,./2_bowtie/seq_hrm/CK_4.hrm.2.fastq,./2_bowtie/seq_hrm/CK_5.hrm.2.fastq,./2_bowtie/seq_hrm/CK_6.hrm.2.fastq,./2_bowtie/seq_hrm/T1_1.hrm.2.fastq,./2_bowtie/seq_hrm/T1_2.hrm.2.fastq,./2_bowtie/seq_hrm/T1_3.hrm.2.fastq,./2_bowtie/seq_hrm/T1_4.hrm.2.fastq,./2_bowtie/seq_hrm/T1_5.hrm.2.fastq,./2_bowtie/seq_hrm/T1_6.hrm.2.fastq,./2_bowtie/seq_hrm/T3_1.hrm.2.fastq,./2_bowtie/seq_hrm/T3_2.hrm.2.fastq,./2_bowtie/seq_hrm/T3_3.hrm.2.fastq,./2_bowtie/seq_hrm/T3_4.hrm.2.fastq,./2_bowtie/seq_hrm/T3_5.hrm.2.fastq,./2_bowtie/seq_hrm/T3_6.hrm.2.fastq,./2_bowtie/seq_hrm/T5_1.hrm.2.fastq,./2_bowtie/seq_hrm/T5_2.hrm.2.fastq,./2_bowtie/seq_hrm/T5_3.hrm.2.fastq,./2_bowtie/seq_hrm/T5_4.hrm.2.fastq,./2_bowtie/seq_hrm/T5_5.hrm.2.fastq,./2_bowtie/seq_hrm/T5_6.hrm.2.fastq,./2_bowtie/seq_hrm/T6_1.hrm.2.fastq,./2_bowtie/seq_hrm/T6_2.hrm.2.fastq,./2_bowtie/seq_hrm/T6_3.hrm.2.fastq,./2_bowtie/seq_hrm/T6_4.hrm.2.fastq,./2_bowtie/seq_hrm/T6_5.hrm.2.fastq,./2_bowtie/seq_hrm/T6_6.hrm.2.fastq \
  --presets meta-sensitive \
  --memory 0.9 \
  --mem-flag 1 \
  --num-cpu-threads 72 \
  --out-dir ./6_assembly/2_megahit_plus \
  --out-prefix mg_6 \
  --keep-tmp-files \
  --tmp-dir /database_new/work/zryan/TEMP_analysis/mg_6_megahit

# 验证结果
if [ ! -f ./6_assembly/2_megahit_plus/mg_6_plus.contigs.fa ]; then
  echo "❌ MEGAHIT 组装失败！"
  exit 1
fi

echo "✅ MEGAHIT 完成！结果: ./6_assembly/2_megahit_plus/mg_6_plus.contigs.fa"



## 子任务5：宏基因组分箱binning（使用MetaBAT2生成bins）
# 输出很多bin.数字.fa文件，每一个fasta相当于一个物种草图


#!/bin/bash
# 脚本名称：11_metabat2_binning_custom.sh
# 功能：基于统一 contigs 和各样本去宿主 reads 进行 MetaBAT2 分箱
# 输入：
#   - Contigs: ./D_3_contigs/mg_6.contigs.fa
#   - 去宿主 reads: ./D_2_HRreads/cleanreads/{sample}.hrm.1/2.fastq
# 输出：
#   - MAGs: ./D_4_bins/bin.*

set -e

# 路径定义
contigs_fa="./D_3_contigs/mg_6.contigs.fa"
reads_dir="./D_2_HRreads/cleanreads"
output_dir="./D_4_bins"
tmp_dir="./D_4_bins/tmp_dir"

# 创建输出与临时目录
mkdir -p "${output_dir}" "${tmp_dir}"

# 1. 构建 Bowtie2 索引（基于你的 contigs）
echo "🔍 正在构建 Bowtie2 索引..."
bowtie2-build "${contigs_fa}" "${tmp_dir}/contigs"

# 2. 定义单样本比对函数
process_sample() {
  local sample_id=$1
  local reads_dir=$2
  local tmp_dir=$3

  # 输入文件
  local r1="${reads_dir}/${sample_id}.hrm.1.fastq"
  local r2="${reads_dir}/${sample_id}.hrm.2.fastq"
  local sam_out="${tmp_dir}/${sample_id}.sam"
  local bam_out="${tmp_dir}/${sample_id}.bam"

  # 比对
  bowtie2 -x "${tmp_dir}/contigs" -1 "${r1}" -2 "${r2}" \
          -p 12 --very-sensitive -S "${sam_out}"

  # 转 BAM + 排序
  samtools view -@ 6 -bS "${sam_out}" | samtools sort -@ 6 -o "${bam_out}"

  # 删除 SAM 以节省空间
  rm -f "${sam_out}"
}

# 导出函数供 parallel 使用
export -f process_sample

# 3. 自动提取样本 ID 列表（从 .hrm.1.fastq 文件）
samples=$(ls "${reads_dir}"/*.hrm.1.fastq | xargs -n1 basename | sed 's/\.hrm\.1\.fastq$//')

# 4. 并行比对所有样本（使用 6 个并行任务，每个用 12 线程，总线程 ≤72）
echo "🧬 正在并行比对所有样本到 contigs..."
parallel -j 6 process_sample {} "${reads_dir}" "${tmp_dir}" ::: $samples

# 5. 合并 BAM 覆盖度生成 depth.txt
echo "📊 正在汇总覆盖度..."
jgi_summarize_bam_contig_depths --outputDepth "${tmp_dir}/depth.txt" "${tmp_dir}"/*.bam

# 6. 运行 MetaBAT2 分箱
echo "📦 正在运行 MetaBAT2 分箱..."
metabat2 -t 72 -i "${contigs_fa}" -a "${tmp_dir}/depth.txt" -o "${output_dir}/bin"

# 7. 清理临时 BAM 文件（可选，若磁盘紧张）
# rm -f "${tmp_dir}"/*.bam

echo "✅ MetaBAT2 分箱完成！"
echo "📁 MAGs 位于: ${output_dir}/"
echo "   文件名示例: bin.0.fa, bin.1.fa, ..."


## 子任务6：使用prodigal得到ORFs（CDS）


#!/bin/bash
# 脚本名称：6_prodigal_orf.sh
# 功能：从组装 contigs 预测 ORFs（CDS），用于后续功能注释
# 前置数据准备：
#   - contigs.fa（来自 MEGAHIT 或 SPAdes）
# 输出：
#   - output_dir/genes.gff  # 包含ORF的 基因组坐标信息（起始位点、终止位点、链方向等）
#   - output_dir/genes.fasta  # 包含ORF的核苷酸序列（即 DNA 序列）
#   - output_dir/proteins.faa  # 包含ORF的蛋白质序列（即氨基酸序列）

set -e

output_dir="output_dir"
mkdir -p "${output_dir}"

prodigal -i contigs.fa \
         -o "${output_dir}/genes.gff" \
         -d "${output_dir}/genes.fasta" \
         -a "${output_dir}/proteins.faa" \
         -p meta -f gff

echo "✅ Prodigal ORF 预测完成"



## 子任务7：构建非冗余基因目录（用于丰度定量）
#!/bin/bash
# 脚本名称：10_cdhit_nonredundant.sh
# 功能：构建非冗余基因目录（用于丰度定量）
# 前置数据准备：
#   - genes.fasta（来自 Prodigal）
# 输出：
#   - output_dir/nonredundant_genes.fasta

# cdhit参数
# -i：fasta格式的输入序列文件，多个宏基因组的基因序列需要合并到一起
# -o：输出文件的文件名
# -c：序列相似度identity阈值，默认为0.9
# -G：设置全局比对还是局部比对，默认为1也即全局比对，如果设置0也即局部比对，需要配和覆盖率参数使-A、-aL、AL、-aS、-AS、-U、-uL、-uS
# -M：内存限制(MB)，默认为800，设置0则无限制
# -T：程序运行使用的核数，默认为1，设置为0则使用所有核
# -n：word过滤时的word长度，默认为5，具体如下所示：
# -n 5 for -c 0.7 ~ 1.0
# -n 4 for -c 0.6 ~ 0.7
# -n 3 for -c 0.5 ~ 0.6
# -n 2 for -c 0.4 ~ 0.5
# -l：分析序列的最短长度，低于此长度的序列被丢掉，默认为10
# -t：对于冗余的容忍度，默认为2，也即去冗余后还可能会保留有2%的冗余
# -d：聚类信息文件中各个聚类组中序列名的长度，默认为20，设为0则将取完整序列名
# -s：序列长度差异阈值，默认为0，如果设置0.9较短序列应该达到代表序列长度的90%
# -S：序列长度差异阈值，默认为999999，如果设置为60，较短序列与代表序列的长度差异不能超过69个氨基酸
# -aL：控制代表序列比对覆盖率的参数，默认为0，若设为0.9则表示比对区间要占到较长序列的90%
# -AL：控制代表序列比对覆盖率的参数，默认为99999999，如果设置为60，比对的序列中较长序列长度为400，那么比对长度应大于340
# -aS：控制代表序列比对覆盖率的参数，默认为0，如果设置为0.9，那么比对区间应占到较短序列长度的90%
# -AS：控制代表序列比对覆盖率的参数，默认为99999999，如果设置为60，比对的序列中较短序列长度为400，那么比对长度应大于340
# -A：两条序列最小的比对覆盖率，默认为0
# -uL：对较长序列最大不匹配的比例，默认为1.0，如果设置为0.1，不匹配区间不能超过较长序列的10%
# -uS：对较短序列最大不匹配的比例，默认为1.0，如果设置为0.1，不匹配区间不能超过较短序列的10%
# -U：最大的不匹配长度，默认为99999999
# -p：默认为0，设置为1则在聚类文件中打印详细的比对情况
# -g：是否开启精确模式，默认为0也即关闭。在默认算法中，一个序列会依次与代表序列进行比对直到满足相似度阈值，而设置为1则会与所有代表序列进行比对，选择最佳的相似度进行聚类
# -sc：默认为0，也即根据代表序列长度对聚类簇进行排序，设置为1则根据聚类簇的大小（也即每个聚类簇的序列数目）进行排序
# -sf：默认为0，也即根据代表序列长度对输出fasta序列，设置为1则根据聚类簇的大小（也即每个聚类簇的序列数目）对输出序列进行排序

set -e

output_dir="output_dir"
mkdir -p "${output_dir}"

cd-hit-est \
  -i genes.fasta \
  -o "${output_dir}/nonredundant_genes.fasta" \
  -M 0 -T 16

echo "✅ 非冗余基因集构建完成"



## 子任务8：将MAGs转化为MAG水平蛋白质序列
# 目的：为后续功能注释（如 ARGs、MGEs、KEGG、COG 等）提供蛋白质序列输入。

#!/bin/bash
# 脚本名称：18_prodigal_mag_to_proteins.sh
# 功能：对 mag_directory/ 中的每个 MAG 运行 Prodigal，预测蛋白质序列
# 前置数据准备：
#   - mag_directory/ （包含 *.fa 的 MAGs）
# 输出：
#   - mag_proteins/ （每个 MAG 对应一个 .faa 蛋白质文件）

set -e

mag_directory="./7_binning"
mag_proteins_dir="./10_mags_protein"
tmp_dir="./10_mags_protein/tmp"

# 创建输出目录
mkdir -p "${mag_proteins_dir}"
mkdir -p "${tmp_dir}"

# 导出输出目录，使函数内可访问
export mag_proteins_dir

# 定义单个 MAG 的蛋白质预测函数
predict_protein() {
  local mag_file=$1
  local base_name=$(basename "$mag_file" .fa)
  local protein_out="${mag_proteins_dir}/${base_name}.faa"

  prodigal \
    -i "$mag_file" \
    -a "$protein_out" \
    -p meta \
    -q
}

# 导出函数供 parallel 使用
export -f predict_protein

# 获取所有 MAG 文件列表（安全方式：使用数组或直接通配）
# 使用 printf 避免 ls 在极端情况下的问题（虽此处影响不大）
mag_files=( "${mag_directory}"/*.fa )

# 并行运行 Prodigal
parallel -j 24 predict_protein ::: "${mag_files[@]}"

echo "✅ 所有 MAG 蛋白质预测完成！结果在 ${mag_proteins_dir}/ (*.faa)"



## 子任务9、过滤低质量 MAGs（checkM2）
# 分析流程：在获得 bins 之后，功能注释或 GTDB-Tk 分析之前。
# 目的：根据 完整性（completeness） 和 污染度（contamination） 筛选得到高质量 MAGs，通常标准为：Completeness ≥ 70%；Contamination ≤ 10%
# 后续用途：将 filtered_mags/ 作为 GTDB-Tk 注释、ARG/MGE 分析、丰度计算 的新输入目录，替代原始 mag_directory/，以确保下游分析基于可靠基因组。

## checkm2环境的准备
# 注意:checkm2不能直接mamba install,会解析错误,需要用git clone的方式进行
# 下载配置文件,编译所需文件
git clone --recursive https://github.com/chklovski/checkm2.git && cd checkm2
# 创建conda环境,并安装所有checkm2的依赖
mamba env create -n checkm2 -f checkm2.yml -y
# 激活环境
conda activate checkm2
# 安装checkm2本体
python setup.py install
# 检查是否正常安装(输出帮助)
checkm2 -h
# 安装数据库(直接安装在自定义文件夹,自己标注版本号)
checkm2 database --download --path /database/work/zryan/biodb/checkm2_1.1.0


## checkm2脚本

#!/bin/bash
# 脚本名称：19_checkm2_filter_mags.sh
# 功能：使用 CheckM2 评估 MAG 质量，并筛选高质量 MAGs（Completeness ≥70%, Contamination ≤10%）
# 输入：
#   - 原始 MAGs 目录：/database/work/zryan/analysis/mg_6_new/D_5_MAGs/
#   - MAG 文件命名格式：bin.*.fa （如 bin.1.fa, bin.23.fa）
# 输出（统一在 ./D_9_filtered_MAGs/ 下）：
#   - checkm_output/        → CheckM2 原始输出（含 quality_report.tsv）
#   - filtered_mags/        → 符合质量标准的 MAGs (.fa)
#   - checkm_filtered_quality.tsv → 筛选后的质量报告

set -e

# === 配置区 ===
INPUT_MAG_DIR="./D_4_bins"
OUTPUT_BASE="./D_5_MAGs"
# 数据库：需要定位到文件而不是所在路径
CHECKM2_DB_PATH="/database/work/zryan/biodb/checkm2_1.1.0/CheckM2_database/uniref100.KO.1.dmnd"

# 创建输出目录
mkdir -p "${OUTPUT_BASE}/checkm_output" "${OUTPUT_BASE}/filtered_mags"

# === Step 1: 运行 CheckM2 质量预测 ===
echo "🔬 正在运行 CheckM2 质量评估..."
checkm2 predict \
  --threads 72 \
  --input "${INPUT_MAG_DIR}" \
  --output-directory "${OUTPUT_BASE}/checkm_output" \
  --database_path "${CHECKM2_DB_PATH}" \
  --extension ".fa" \
  --force

# === Step 2: 筛选高质量 MAGs ===
QUALITY_REPORT="${OUTPUT_BASE}/checkm_output/quality_report.tsv"
FILTERED_TSV="${OUTPUT_BASE}/checkm_filtered_quality.tsv"

# 列顺序（根据 CheckM2 v1.1.0 输出）：
# Column 1: Name
# Column 2: Completeness
# Column 3: Contamination
# 所以条件：$2 >= 50 && $3 <= 10
awk -F'\t' 'NR==1 || ($2 >= 50 && $3 <= 10)' "${QUALITY_REPORT}" > "${FILTERED_TSV}"

# === Step 3: 复制合格 MAGs 到 filtered_mags/ ===
echo "📦 正在复制高质量 MAGs..."
while IFS= read -r bin_id; do
  if [ -n "$bin_id" ] && [ "$bin_id" != "Name" ]; then
    SRC_FILE="${INPUT_MAG_DIR}/${bin_id}.fa"
    if [ -f "$SRC_FILE" ]; then
      cp "$SRC_FILE" "${OUTPUT_BASE}/filtered_mags/"
    else
      echo "⚠️ 警告：未找到 ${SRC_FILE}，跳过。"
    fi
  fi
done < <(tail -n +2 "${FILTERED_TSV}" | cut -f1)

# === 完成提示 ===
TOTAL_MGS=$(wc -l < "${QUALITY_REPORT}" | awk '{print $1-1}')
FILTERED_COUNT=$(wc -l < "${FILTERED_TSV}" | awk '{print $1-1}')

echo ""
echo "✅ CheckM2 质控完成！"
echo "  - 总 MAGs 数量: ${TOTAL_MGS}"
echo "  - 高质量 MAGs (C≥70%, Cont≤10%): ${FILTERED_COUNT}"
echo "  - 质量报告: ${FILTERED_TSV}"
echo "  - 高质量 MAGs 已保存至: ${OUTPUT_BASE}/filtered_mags/"



## 子任务30：使用dRep工具，将MAGs进行去冗余，得到物种水平的代表性基因组（SRGs）


#!/bin/bash
# 脚本名称：4_drep_srg.sh
# 功能：对 filtered MAGs 进行去冗余，生成 SRGs
# 输入：
#   - ./4_filtered_MAGs/fbin.*.fa
# 输出目录：./1_drep_out/
# 输出文件：./1_drep_out/dereplicated_genomes/*.fna （即 SRGs）
# 笔记:
# （1）目的：对多个样本的MAGs进行去冗余，生成非冗余的物种代表基因组（SRGs）。是完整流程（聚类+选择最优基因组）。
# （2）适用性：广泛用于宏基因组MAGs后处理，适用于肠道等多样本研究。
# （3）归类：第一部分：数据准备模块（按题设）
# （4）工具与环境：drep，conda可安装，Linux支持良好。
# （5）数据准备：filtered MAGs（fbin.1.fa, fbin.2.fa...）



set -e

INPUT_DIR="./4_filtered_MAGs"
OUTPUT_DIR="./1_drep_out"

dRep dereplicate "$OUTPUT_DIR" \
    -g "$INPUT_DIR"/fbin.*.fa \
    -p 32 \
    --completeness 50 \
    --contamination 10 \
    --S_algorithm fastANI \
    --S_ani 0.95 \
    --cov_thresh 0.9

echo "✅ dRep 完成！SRGs 保存至 $OUTPUT_DIR/dereplicated_genomes/"






# 第二部分：基于reads层面的分析
## 子任务29：使用AGS-and-ACN-tools计算平均基因组大小（AGS）和16S rRNA基因平均拷贝数（ACN）

#!/bin/bash
# 脚本名称：3_ags_acn.sh
# 功能：计算 AGS 和 ACN
# 输入：
#   - ./1_kneaddata/sample1_hrm_1.fastq, ./1_kneaddata/sample1_hrm_2.fastq ...
# 工具路径：
#   - AGS_and_ACN_dir=/tools/AGS-and-ACN-tools
# 参考16S数据库（内置）：
#   - SILVA_16S.fasta 在工具包内
# 输出目录：./2_ags_acn/
# 笔记:
# （1）目的：基于未组装的宏基因组reads估算群落水平的平均基因组大小与16S拷贝数，用于生态推断。是单步骤分析。
# （2）适用性：专为宏基因组设计，适用于所有微生物组场景（包括肠道），且必须使用 raw 或 clean reads（未组装）。
# （3）归类：第二部分：基于reads层面的分析
# （4）工具与环境：AGS-and-ACN-tools（Python脚本），依赖 bbmap 和 bwa。可构建 conda 环境运行。
# （5）数据准备：clean reads（sample1.hrm.1.fastq, sample1.hrm.2.fastq...）

set -e

AGS_DIR="/tools/AGS-and-ACN-tools"
OUT_DIR="./2_ags_acn"
mkdir -p "$OUT_DIR"

# 获取样本列表
samples=$(ls ./1_kneaddata/*_hrm_1.fastq | xargs -n1 basename | sed 's/_hrm_1.fastq//')

for sample in $samples; do
    echo "Running AGS/ACN for $sample ..."
    python "$AGS_DIR"/run_AGS_ACN.py \
        --forward ./1_kneaddata/"${sample}"_hrm_1.fastq \
        --reverse ./1_kneaddata/"${sample}"_hrm_2.fastq \
        --output "$OUT_DIR"/"${sample}"_ags_acn.txt \
        --threads 8
done

echo "✅ AGS/ACN 计算完成！结果保存至 $OUT_DIR"




## 子任务10：使用kraken和bracken进行微生物分类和丰度定量

#!/bin/bash
# 脚本名称：3_kraken2_bracken.sh
# 功能：微生物群落物种组成分析
# 前置数据准备：
#   - 去宿主 reads: sample1.hrm.1.fastq ... sample4.hrm.2.fastq
#   - Kraken2 数据库路径已知（此处假设为 KRAKEN_DB）
# 输出：
#   - output_dir/*.kreport
#   - output_dir/*.S.bracken
#   - output_dir/*.S.bracken.kreport

set -e

output_dir="output_dir"
mkdir -p "${output_dir}"
KRAKEN_DB="/path/to/kraken2_database"  # ← 用户需替换

process_sample() {
  local sample=$1
  kraken2 --threads 8 --db "$KRAKEN_DB" \
          --paired "${sample}.hrm.1.fastq" "${sample}.hrm.2.fastq" \
          --report "${output_dir}/${sample}.kreport" \
          --memory-mapping

  bracken -d "$KRAKEN_DB" -i "${output_dir}/${sample}.kreport" \
          -o "${output_dir}/${sample}.S.bracken" \
          -w "${output_dir}/${sample}.S.bracken.kreport" \
          -l S -t 8
}

export -f process_sample

samples=("sample1" "sample2" "sample3" "sample4")
parallel -j 4 process_sample ::: "${samples[@]}"

echo "✅ Kraken2 + Bracken 分析完成"


## 子任务11：ARGs_OAP 抗性基因注释（reads层面）

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

## 子任务12：HUMAnN3功能注释

### 多样本脚本（正常用这个）####

#!/bin/bash
# 脚本名称：5_humann3_reads.sh
# 功能：基于去宿主 reads 进行通路级功能注释
# 前置数据准备：
#   - 去宿主 reads: sample1.hrm.1.fastq ... sample4.hrm.2.fastq
#   - HUMAnN3 及 MetaPhlAn 数据库已安装
# 输出：
#   - output_dir/{genefamilies,pathabundance} 归一化表格

# HUMAnN3 批处理脚本（20251130）—— 修复 MetaPhlAn 数据库兼容性问题
# 功能：基于去宿主后的双端 fastq，自动合并为单端，调用 HUMAnN3 进行功能注释，并归一化
# 并行数：64（预留8核给系统）
# 输出目录：./5_humann3/
# 临时目录：/database_new/work/zryan/TEMP_analysis/mg_6_humann3

set -e  # 遇错退出

# 1. 创建输出与临时目录
mkdir -p ./5_humann3/{gene_families,path_abundance,path_coverage,reports}
mkdir -p /database_new/work/zryan/TEMP_analysis/mg_6_humann3

# 2. 定义单样本处理函数
process_sample() {
    local sample=$1
    local tmpdir="/database_new/work/zryan/TEMP_analysis/mg_6_humann3"
    local outdir="./5_humann3"
 
    # 输入文件路径（注意：原始文件是 .hrm.1.fastq / .hrm.2.fastq）
    local r1_in="./2_bowtie/seq_hrm/${sample}.hrm.1.fastq"
    local r2_in="./2_bowtie/seq_hrm/${sample}.hrm.2.fastq"
    
    # 临时单端文件（直接写入 tmpdir，不复制原始文件）
    local merged_fastq="${tmpdir}/${sample}.fastq"

    # 检查输入是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "⚠️  Warning: Missing input for sample $sample. Skipping."
        return 0
    fi

    # 直接合并双端为单端（无需先 cp 到 tmpdir）
    cat "$r1_in" "$r2_in" > "$merged_fastq"

    # HUMAnN3 输出子目录（每个样本独立）(这一步产生的数据会产生大量样本.)
    local sample_out_dir="${outdir}/reports/${sample}_humann3_out"
    mkdir -p "$sample_out_dir"

    # === 关键修改：显式导出 MetaPhlAn 数据库环境变量 ===
    export METAPHLAN_DB_DIR="/database/work/zryan/conda_envs/humann3/lib/python3.8/site-packages/metaphlan/metaphlan_databases"
    export METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202307"

    # 运行 HUMAnN3（使用合并后的单端文件）
    humann \
        --input "$merged_fastq" \
        --output "$sample_out_dir" \
        --threads 12 \
        --memory-use maximum \
        --search-mode uniref90 \
        --remove-temp-output

    # === Gene Families ===
    humann_renorm_table \
        --input "${sample_out_dir}/${sample}_genefamilies.tsv" \
        --output "${outdir}/gene_families/${sample}.hrm_genefamilies_relab.tsv" \
        --units relab

    humann_renorm_table \
        --input "${sample_out_dir}/${sample}_genefamilies.tsv" \
        --output "${outdir}/gene_families/${sample}.hrm_genefamilies_cpm.tsv" \
        --units cpm

    # === Pathway Abundance ===
    humann_renorm_table \
        --input "${sample_out_dir}/${sample}_pathabundance.tsv" \
        --output "${outdir}/path_abundance/${sample}.hrm_pathabundance_relab.tsv" \
        --units relab

    humann_renorm_table \
        --input "${sample_out_dir}/${sample}_pathabundance.tsv" \
        --output "${outdir}/path_abundance/${sample}.hrm_pathabundance_cpm.tsv" \
        --units cpm

    # === Pathway Coverage ===（不归一化，直接复制）
    cp "${sample_out_dir}/${sample}_pathcoverage.tsv" \
       "${outdir}/path_coverage/${sample}.hrm_pathcoverage_final.tsv"

    # 清理临时合并文件
    rm -f "$merged_fastq"

    echo "✅ Finished sample: $sample"
}

# 3. 导出函数供 parallel 使用
export -f process_sample

# 4. 提取样本名列表（从 .hrm.1.fastq 文件推断）
samples=$(ls ./2_bowtie/seq_hrm/*.hrm.1.fastq 2>/dev/null | xargs -n1 basename | sed 's/\.hrm\.1\.fastq$//')

if [[ -z "$samples" ]]; then
    echo "❌ No input files found in ./2_bowtie/seq_hrm/ matching pattern *.hrm.1.fastq"
    exit 1
fi

# 5. 并行处理（使用6个任务并行，每个任务用12线程，总计72核以内）
echo "🚀 Starting HUMAnN3 on $(echo $samples | wc -w) samples using 6 parallel jobs (12 threads each)..."
parallel -j 6 process_sample ::: $samples

# 6. 汇总所有样本的表格（按类别合并）
echo "📊 Merging tables..."

# 基因家族（relab）
humann_join_tables -i ./5_humann3/gene_families/*_genefamilies_relab.tsv -o ./5_humann3/gene_families_relab_merged.tsv
humann_join_tables -i ./5_humann3/gene_families/*_genefamilies_cpm.tsv -o ./5_humann3/gene_families_cpm_merged.tsv

# 通路丰度（relab）
humann_join_tables -i ./5_humann3/path_abundance/*_pathabundance_relab.tsv -o ./5_humann3/path_abundance_relab_merged.tsv
humann_join_tables -i ./5_humann3/path_abundance/*_pathabundance_cpm.tsv -o ./5_humann3/path_abundance_cpm_merged.tsv

# 通路覆盖度
humann_join_tables -i ./5_humann3/path_coverage/*_pathcoverage_final.tsv -o ./5_humann3/path_coverage_merged.tsv

# 7. 清理临时目录（可选）
rm -rf /database_new/work/zryan/TEMP_analysis/mg_6_humann3

echo "🎉 All done! Results in ./5_humann3/"

### 单样本脚本（如果出现报错，单独运行这个把报错的重跑完）####

#!/bin/bash
# HUMAnN3 单样本补充运行脚本（20251206）
# 用途：重跑因 OOM/中断而失败的单个样本
# 输入：一个样本名（如 CK_2）,如果需要开多个样本,则开多个session就可以.
# 用法: bash ./5_humann3_singlesample.sh CK_2
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_name>"
    echo "Example: $0 CK_2"
    exit 1
fi
sample=$1
echo "🔧 Rescuing sample: $sample"
# === 配置路径（与主脚本一致）===
tmpdir="/database_new/work/zryan/TEMP_analysis/mg_6_humann3"
outdir="./5_humann3"
r1_in="./2_bowtie/seq_hrm/${sample}.hrm.1.fastq"
r2_in="./2_bowtie/seq_hrm/${sample}.hrm.2.fastq"
merged_fastq="${tmpdir}/${sample}.fastq"
sample_out_dir="${outdir}/reports/${sample}_humann3_out"

# === 检查输入是否存在 ===
if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
    echo "❌ Error: Input files not found for sample $sample. Exiting."
    exit 1
fi

# === 1. 清理旧的不完整中间目录 ===
if [ -d "$sample_out_dir" ]; then
    echo "🗑️  Removing existing incomplete output directory: $sample_out_dir"
    rm -rf "$sample_out_dir"
fi

# === 2. 创建所需目录 ===
mkdir -p "$outdir"/{gene_families,path_abundance,path_coverage,reports}
mkdir -p "$tmpdir"

# === 3. 合并双端为单端（与你当前流程一致）===
echo "🔀 Merging R1 and R2 into single-end FASTQ..."
cat "$r1_in" "$r2_in" > "$merged_fastq"

# === 4. 设置 MetaPhlAn 数据库环境变量 ===
export METAPHLAN_DB_DIR="/database/work/zryan/conda_envs/humann3/lib/python3.8/site-packages/metaphlan/metaphlan_databases"
export METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202307"

# === 5. 运行 HUMAnN3 ===
echo "🚀 Running HUMAnN3 for $sample..."
humann \
    --input "$merged_fastq" \
    --output "$sample_out_dir" \
    --threads 72 \
    --memory-use maximum \
    --search-mode uniref90 \
    --remove-temp-output

# === 6. 归一化与复制 ===
echo "📊 Normalizing and copying results..."

humann_renorm_table \
    --input "${sample_out_dir}/${sample}_genefamilies.tsv" \
    --output "${outdir}/gene_families/${sample}.hrm_genefamilies_relab.tsv" \
    --units relab

humann_renorm_table \
    --input "${sample_out_dir}/${sample}_genefamilies.tsv" \
    --output "${outdir}/gene_families/${sample}.hrm_genefamilies_cpm.tsv" \
    --units cpm

humann_renorm_table \
    --input "${sample_out_dir}/${sample}_pathabundance.tsv" \
    --output "${outdir}/path_abundance/${sample}.hrm_pathabundance_relab.tsv" \
    --units relab

humann_renorm_table \
    --input "${sample_out_dir}/${sample}_pathabundance.tsv" \
    --output "${outdir}/path_abundance/${sample}.hrm_pathabundance_cpm.tsv" \
    --units cpm

cp "${sample_out_dir}/${sample}_pathcoverage.tsv" \
   "${outdir}/path_coverage/${sample}.hrm_pathcoverage_final.tsv"

# === 7. 清理临时合并文件 ===
rm -f "$merged_fastq"

echo "✅ Successfully rescued sample: $sample"
echo "📁 Final outputs in ./5_humann3/{gene_families,path_abundance,path_coverage}/"





第三部分：基于contigs、ORF层面的分析
## 子任务45：使用Quast工具，对基因组组装结果进行质量评估
#!/bin/bash
# 脚本名称：19_quast_contigs.sh
# 功能：组装质量评估
# 输入：
#   - ./2_assembly/contigs.fa
# 输出目录：./3_quast_out/
# 笔记：
# （1）目的：评估 contigs 或 MAGs 的组装质量（N50, #contigs, completeness 等）。
# （2）适用性：通用。
# （3）归类：第三部分（若输入 contigs）或第四部分（若输入 MAGs）→ 此处以 contigs 为例 → 第三部分
# （4）工具与环境：quast，conda 安装。
# （5）数据准备：contigs.fa

set -e

quast.py ./2_assembly/contigs.fa -o ./3_quast_out -t 16

echo "✅ QUAST 完成！报告在 ./3_quast_out/report.html"


## 子任务13：reads水平基因丰度定量（样本的回帖）

# 基于非冗余基因集+样本clean reads得到。
# 这一步非常关键：这是所有基于contigs个性化分析最终回归到样本丰度的中间步骤。因为contigs，ORFs等最终只会得到一个全样本的注释结果，必须有回帖的数据才能够得到样本的表达量。
# 注意：本步骤和前面的cd-hit高度绑定。一旦生成 nr_genes.ffn 和 nr_genes.faa，其序列和 ID 不能再改动（不能重新 CD-HIT、不能删减）

### 版本1:准比对法（Pseudo-alignment）（salmon，通用）

#!/bin/bash
# 脚本名称：12_salmon_quant_genes.sh
# 功能：使用 Salmon 非比对方法定量非冗余基因集在各样本中的丰度（TPM + counts）
# 前置数据准备：
#   - nonredundant_genes.fasta（非冗余基因核苷酸序列）
#   - ./clean_reads/ 目录下各样本的 clean reads: sample1.R1.qualified.fastq.gz, sample1.R2.qualified.fastq.gz 等
# 输出：
#   - salmon_index/
#   - quant_results/ （每个样本一个子目录）
#   - gene_counts.txt, gene_tpm.txt（合并后的丰度矩阵）

set -e

# 创建输出目录
mkdir -p salmon_index quant_results merged_output

# 1. 构建 Salmon 索引
salmon index -t nonredundant_genes.fasta -i salmon_index --type quasi -k 31

# 2. 定义处理单个样本的函数
quant_sample() {
  local sample=$1
  local r1="./clean_reads/${sample}.R1.qualified.fastq.gz"
  local r2="./clean_reads/${sample}.R2.qualified.fastq.gz"
  local out_dir="quant_results/${sample}"

  salmon quant -i salmon_index \
               -l A \
               -1 "$r1" -2 "$r2" \
               -p 12 \
               --meta \
               -o "$out_dir"
}

# 导出函数供 parallel 使用
export -f quant_sample

# 获取样本名列表（适配 .gz 后缀）
samples=$(ls ./clean_reads/*.R1.qualified.fastq.gz | xargs -n1 basename | sed 's/\.R1\.qualified\.fastq\.gz//')

# 并行定量（假设服务器有72核，此处用 -j 6，每样本12线程）
parallel -j 6 quant_sample ::: $samples

# 3. 合并所有样本的 counts 和 TPM
salmon quantmerge --quants quant_results/*/ --column numreads -o merged_output/gene_counts.txt
salmon quantmerge --quants quant_results/*/ --column tpm -o merged_output/gene_tpm.txt

echo "✅ Salmon 基因丰度定量完成！结果在 merged_output/"



### 版本2:全比对法（Full alignment）（bwa）
# 注意：此方法计算的是“比对到基因的 reads 数”，但不直接输出 TPM/FPKM，需后续用 bedtools coverage 或自定义脚本计算长度归一化丰度。此处仅完成比对与基础统计。

#!/bin/bash
# 脚本名称：13_bwa_quant_genes.sh
# 功能：使用 BWA 比对 reads 到非冗余基因集，生成 BAM 文件用于丰度估计
# 前置数据准备：
#   - nonredundant_genes.fasta
#   - ./clean_reads/ 下的 clean reads (.fastq.gz)
# 输出：
#   - bwa_index.*
#   - bwa_bams/ （各样本排序后的 BAM）
#   - 注：丰度统计需额外步骤（如 featureCounts 或自定义脚本）

set -e

mkdir -p bwa_index bwa_bams tmp_sam

# 1. 构建 BWA 索引
bwa index -p bwa_index nonredundant_genes.fasta

# 2. 定义单样本处理函数
align_sample() {
  local sample=$1
  local r1="./clean_reads/${sample}.R1.qualified.fastq.gz"
  local r2="./clean_reads/${sample}.R2.qualified.fastq.gz"
  local sam="tmp_sam/${sample}.sam"
  local bam="bwa_bams/${sample}.bam"
  local sorted_bam="bwa_bams/${sample}.sorted.bam"

  # 比对（BWA 支持 .gz 输入）
  bwa mem -t 8 bwa_index "$r1" "$r2" -o "$sam"

  # 转 BAM 并排序
  samtools view -@ 4 -S -b "$sam" | samtools sort -@ 8 -o "$sorted_bam"
  
  # 清理临时 SAM
  rm "$sam"
}

export -f align_sample

samples=$(ls ./clean_reads/*.R1.qualified.fastq.gz | xargs -n1 basename | sed 's/\.R1\.qualified\.fastq\.gz//')

# 并行比对（-j 9，每样本8+4+8≈20线程，总核数控制在72内）
parallel -j 9 align_sample ::: $samples

echo "✅ BWA 比对完成！排序 BAM 在 bwa_bams/。后续可用 bedtools 或 featureCounts 计算丰度。"
echo "💡 建议：若需生成 count 矩阵，推荐使用 featureCounts（来自 Subread 包）而非 bedtools bamtobed，因其可直接按基因 ID 统计 reads 数。"



### 版本3：基于 CoverM 的宏基因组基因丰度定量（宏基因组推荐）

#!/bin/bash
# 脚本名称：14_coverm_quant_genes.sh
# 功能：使用 CoverM 全面定量非冗余基因集在各样本中的多种丰度指标
# 前置数据准备：
#   - 参考基因集: ./D_7_nonredundant_geneset/nonredundant_genes.fasta
#   - clean reads: ./D_2_HRreads/cleanreads/ 下的 *.hrm.1.fastq 和 *.hrm.2.fastq（未压缩）
# 输出：
#   - ./C_1_Gene_quant/coverm_output/ 包含多种丰度矩阵

set -e

# 定义路径
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
GENE_FASTA="./D_7_nonredundant_geneset/nonredundant_genes.fasta"
OUTPUT_DIR="./C_1_Gene_quant/coverm_output"
TMPDIR="./C_1_Gene_quant/tmp"

# 创建输出目录与临时目录（确保存在）
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMPDIR}/allmethods"
mkdir -p "${TMPDIR}/metabat"
mkdir -p "${TMPDIR}/coverage_histogram"

# 运行coverm,其他所有的方法
TMPDIR="${TMPDIR}/allmethods" coverm contig \
  --reference "${GENE_FASTA}" \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --methods mean covered_fraction covered_bases variance count length reads_per_base rpkm tpm \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT" \
  --verbose



# 运行coverm,coverage_histogram方法，因为不能与其他一起运行
TMPDIR="${TMPDIR}/coverage_histogram" coverm contig \
  --reference "${GENE_FASTA}" \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --methods coverage_histogram \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT_coverage_histogram" \
  --verbose




# 运行coverm,单独运行metabat方法.因为不能与其他一起运行
TMPDIR="${TMPDIR}/metabat" coverm contig \
  --reference "${GENE_FASTA}" \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --methods metabat \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT_metabat" \
  --verbose



echo "✅ CoverM 基因丰度定量完成！结果在 ${OUTPUT_DIR}/"
echo "💡 输出包含以下指标：mean, coverage_histogram, covered_fraction, covered_bases, variance, length, count, reads_per_base, rpkm, tpm, metabat"

# 说明：参数设置原因
# --min-read-aligned-length 50 \
# 宏基因组 reads 通常为 150 bp（PE），50 bp 是可靠比对的下限。过低会引入噪声；过高会丢失短 contig 的覆盖信息。类似研究常用 50–75 bp（Sczyrba et al., 2017, *Nature Methods*）。
# --min-read-percent-identity 90 \
# 宏基因组中物种多样性高，90–95% identity 是区分近缘菌株的常用阈值。90% 平衡了灵敏度与特异性，避免将 reads 错配到远缘同源序列（Nayfach et al., 2016, *Genome Research*）。
# --min-read-aligned-percent 80 \
# 	要求至少 80% 的 read 被比对上，防止部分比对（partial alignment）导致的假阳性覆盖。这在短 contigs（<1 kb）中尤为重要（Bishara et al., 2018, *Nature Biotechnology*）。

######end



## 子任务14：毒力因子预测
### VFDB本地化部署（手动运行）####

# 独立环境的conda三部曲
# 环境中一定要有diamond
mamba create -n vfdb
conda activate vfdb
mamba install -c bioconda diamond -y

# 使用diamond的makedb命令创建数据库
# 数据库下载网站（核心/完整的蛋白序列），从浏览器进去下载。
http://www.mgc.ac.cn/VFs/download.htm

# 进入到专门存放数据库的文件夹
cd /database/work/zryan/biodb/VFDB
# 直接下载full的protain的数据库
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
# 解压
gunzip VFDB_setB_pro.fas.gz
# 使用DIAMOND创建数据库
diamond makedb --in VFDB_setB_pro.fas -d VFDB_proteins

### 脚本内容####

#!/bin/bash
# 脚本名称：7_vfdb_virulence.sh
# 功能：基于 ORF 蛋白质序列预测毒力因子
# 前置数据准备：
#   - proteins.faa（来自 Prodigal）
#   - VFDB 数据库已下载并构建 DIAMOND 数据库（vfdb.dmnd）
# 输出：
#   - output_dir/vfdb_hits.tsv

set -e

output_dir="output_dir"
mkdir -p "${output_dir}"

diamond blastp \
  --db VFDB_proteins.dmnd \
  --query proteins.faa \
  --out "${output_dir}/vfdb_hits.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --threads 16 \
  --evalue 1e-5 \
  --id 80 \
  --query-cover 80

echo "✅ VFDB 毒力因子预测完成"


## 子任务50：(后续需R)使用GO数据库，创建中间映射；并以另一个独立脚本给出GO中间映射到注释的操作流程

# 笔记：
# （1）目的：
# - 第一步：从 InterProScan 或其他注释工具输出中提取 Gene Ontology（GO）术语，生成“蛋白ID → GO term”映射表（中间映射）。
# - 第二步：将该映射用于后续丰度汇总或功能富集分析。
# 这是两步流程，第一步为数据转换，第二步为注释整合。
# （2）适用性：适用于需要进行 GO 功能富集分析的宏基因组研究，包括肠道微生物组的功能生态学分析。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：
# - 第一步：使用 awk/grep 在 Linux 中处理 InterProScan 输出（无需额外工具）。
# - 第二步：通常在 R 中完成（因涉及矩阵合并与统计），但可提供 bash 汇总框架。
# 两者均可在 Linux 完成，但第二步强烈建议转 R。
# （5）数据准备：
# - 第一步输入：InterProScan 输出（如 ./3_interpro/proteins.tsv）
# - 第二步输入：蛋白丰度表（如 orf_abundance.tsv） + GO 映射表

# 脚本 24-1：创建 GO 中间映射

#!/bin/bash
# 脚本名称：25_1_go_mapping.sh
# 功能：从 InterProScan TSV 输出提取蛋白到 GO term 的映射
# 输入：
#   - ./3_interpro/proteins.tsv （第14列为 GO terms，以","分隔）
# 输出目录：./3_go/
# 输出文件：./3_go/protein2go.tsv （两列：protein_id \t GO_term）

set -e

mkdir -p ./3_go

awk -F'\t' '
    $14 != "" && $14 != "-" {
        n = split($14, go_array, ",");
        for (i = 1; i <= n; i++) {
            if (go_array[i] != "") print $1 "\t" go_array[i]
        }
    }
' ./3_interpro/proteins.tsv | sort -u > ./3_go/protein2go.tsv

echo "✅ GO 中间映射完成！结果保存至 ./3_go/protein2go.tsv"


# 脚本 24-2：GO 注释整合（示例：按 GO term 汇总 ORF 丰度）

# ⚠️ 注：实际丰度汇总需数值计算，bash 不擅长。此处提供逻辑框架，建议后续用 R。
# ⚠️ 下面是R语言脚本

#!/bin/bash
# 脚本名称：25_2_go_summarize.sh
# 功能：将 ORF 丰度表按 GO term 汇总（仅展示逻辑，实际推荐 R 实现）
# 输入：
#   - ./3_go/protein2go.tsv
#   - ./abundance/orf_abundance.tsv （第一列为 protein_id，其余为样本丰度）
# 输出：./3_go/go_abundance_prelim.tsv （未加权求和，仅连接）

set -e

# 此操作在 bash 中效率低且易错，仅作示意
# 真实分析应使用 R：
#   library(dplyr)
#   go_map <- read.delim("protein2go.tsv", header=F)
#   abund  <- read.delim("orf_abundance.tsv", row.names=1)
#   merged <- left_join(go_map, abund, by=c("V1"="row.names"))
#   go_abund <- merged %>% group_by(V2) %>% summarise(across(where(is.numeric), sum))

echo "⚠️ 建议在 R 中完成 GO 丰度汇总。此脚本仅生成连接表供参考。"

join -1 1 -2 1 <(sort ./3_go/protein2go.tsv) <(sort ./abundance/orf_abundance.tsv) > ./3_go/go_abundance_prelim.tsv

echo "✅ 初步连接完成（非最终丰度）！结果在 ./3_go/go_abundance_prelim.tsv"



## 子任务15：KEGG通路注释

### 准备工作
# conda三步曲
mamba create -n kofamscan
conda activate kofamscan
mamba install kofamscan -c bioconda

# 检查依赖，如果没有这三个，则分别安装
# mamba install -c conda-forge ruby
# mamba install -c bioconda hmmer
# mamba install -c conda-forge parallel

# 切换工作目录
cd /database/work/zryan/biodb/kofamscan_db

# 下载KO list和profile HMM文件(ftp和http都可以)
# wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
# wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

# 解压（不要切换工作目录）
gunzip ko_list.gz
tar -xvzf profiles.tar.gz

# 注释脚本

#!/bin/bash
# 脚本名称：8_kegg_annotation.sh
# 功能：KEGG Orthology 注释（使用 KofamScan）
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - ko_list: /database/work/zryan/biodb/kofamscan_db/ko_list
#   - profiles: /database/work/zryan/biodb/kofamscan_db/profiles （含 .hmm 文件的目录）
# 临时目录：./C_3_KEGG/tmp
# 输出目录：./C_3_KEGG/
# 输出文件：./C_3_KEGG/ko_output.txt, ./C_3_KEGG/ko_output.detail.txt, ./C_3_KEGG/ko_output.detail_tsv.txt

set -e

# 定义路径
input_proteins="./D_6_ORFs/proteins.faa"
ko_list="/database/work/zryan/biodb/kofamscan_db/ko_list"
profiles_dir="/database/work/zryan/biodb/kofamscan_db/profiles"
output_dir="./C_3_KEGG"
tmp_dir="${output_dir}/tmp"

# 创建输出和临时目录
mkdir -p "${output_dir}"
mkdir -p "${tmp_dir}"

# 第一次：生成 mapper 格式（用于通路分析）
exec_annotation "${input_proteins}" \
  -f mapper \
  -p "${profiles_dir}" \
  -k "${ko_list}" \
  --cpu 72 \
  --tmp-dir "${tmp_dir}/mapper" \
  -o "${output_dir}/ko_output.mapper.txt"

# 第二次：生成 detail 格式（用于结果审查或自定义过滤）
exec_annotation "${input_proteins}" \
  -f detail \
  -p "${profiles_dir}" \
  -k "${ko_list}" \
  --cpu 72 \
  --tmp-dir "${tmp_dir}/detail" \
  -o "${output_dir}/ko_output.detail.txt"

  # 第三次：生成 detail-tsv 格式（用于结果审查或自定义过滤）
exec_annotation "${input_proteins}" \
  -f detail-tsv \
  -p "${profiles_dir}" \
  -k "${ko_list}" \
  --cpu 72 \
  --tmp-dir "${tmp_dir}/detail_tsv" \
  -o "${output_dir}/ko_output.detail_tsv.txt"

echo "✅ KEGG 注释完成！结果保存至 ${output_dir}/ko_output.txt"


## 子任务47：使用MetaCyC数据库进行代谢通路注释

#!/bin/bash
# 脚本名称：21_metacyc.sh
# 功能：MetaCyc 通路注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - /database/MetaCyc/MetaCyc.dmnd
# 输出：./3_metacyc/metacyc_hits.tsv
# 笔记：
# （1）目的：将基因映射到 MetaCyc 代谢通路。
# （2）适用性：适用于代谢潜力分析。
# （3）归类：第三部分
# （4）工具与环境：通常通过 blastp + 自定义脚本，或使用 Pathway Tools。此处用 DIAMOND。
# （5）数据准备：proteins.faa

set -e

mkdir -p ./3_metacyc

diamond blastp \
    --db /database/MetaCyc/MetaCyc.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_metacyc/metacyc_hits.tsv \
    --outfmt 6 qseqid sseqid pident evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5

echo "✅ MetaCyc 注释完成！结果在 ./3_metacyc/metacyc_hits.tsv"


## 子任务31：使用mobileOG-db工具，将contigs级别的蛋白质序列（ORFs）进行MGEs注释）
# 可以用来补充MAGs水平上的注释。把为分箱的那部分也给补充。


#!/bin/bash
# 脚本名称：5_mobileog.sh
# 功能：使用 mobileOG-db 对 ORFs 进行 MGE 注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库路径：
#   - /database/mobileOG-db/mobileOG.db
# DIAMOND 数据库已建好：/database/mobileOG-db/mobileOG.dmnd
# 输出目录：./3_mobileog/
# 输出文件：mobileog_hits.tsv
# 笔记:
# （1）目的：在ORF层面注释移动遗传元件（MGEs）相关基因。是单步骤功能注释。
# （2）适用性：适用于宏基因组ORFs，包括肠道微生物组。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：mobileOG-db（依赖 DIAMOND + Prodigal），需 Python 3.6.15 环境，可 conda 构建。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_mobileog

diamond blastp \
    --db /database/mobileOG-db/mobileOG.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --tnpdir ./3_mobileog/tnpdir \
    --out ./3_mobileog/mobileog_hits.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 72 \
    --evalue 1e-6 \
    --pidentvalue 80 \
    --queryscore 80 \
    --ultra-sensitive \
    --max-target-seqs 1 \
    --header simple \
    --verbose
  

echo "✅ mobileOG 注释完成！结果保存至 ./3_mobileog/mobileog_hits.tsv"


## 子任务46：使用TCDB数据库对微生物组转运功能基因进行注释

#!/bin/bash
# 脚本名称：20_tcdb.sh
# 功能：转运蛋白注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - /database/TCDB/tcdb.dmnd
# 输出：./3_tcdb/tcdb_hits.tsv
# 笔记：
# （1）目的：注释转运蛋白（Transporter Classification Database）。
# （2）适用性：通用。
# （3）归类：第三部分
# （4）工具与环境：DIAMOND + TCDB，conda 可行。
# （5）数据准备：proteins.faa

set -e

mkdir -p ./3_tcdb

diamond blastp \
    --db /database/TCDB/tcdb.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_tcdb/tcdb_hits.tsv \
    --outfmt 6 qseqid sseqid pident evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5

echo "✅ TCDB 注释完成！结果在 ./3_tcdb/tcdb_hits.tsv"


## 子任务40：使用CAZyme数据库进行碳降解基因注释

#!/bin/bash
# 脚本名称：14_dbcan_cazy.sh
# 功能：CAZyme 注释（使用 dbCAN2 HMM 方法）
# 输入：
#   - ./D_6_ORFs/proteins.faa
# HMM 数据库：
#   - /database/dbCAN2/dbCAN-HMMdb-V11.txt.hmm
# 输出目录：./3_cazy/
# 笔记：
# （1）目的：注释碳水化合物活性酶（CAZymes）。是标准功能注释步骤。
# （2）适用性：适用于所有宏基因组，尤其肠道、土壤等富含纤维素降解菌的环境。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：hmmscan + dbCAN2 HMMs，或 DIAMOND + CAZy。推荐 dbCAN2，conda 安装。
# （5）数据准备：proteins.faa


set -e

mkdir -p ./3_cazy

run_dbcan.py ./D_6_ORFs/proteins.faa protein \
    --out_dir ./3_cazy \
    --db_dir /database/dbCAN2 \
    --dia_cpu 0 --hmm_cpu 32 --hotpep_cpu 0

echo "✅ CAZyme 注释完成！结果在 ./3_cazy/"


## 子任务39：使用NCycDB进行氮循环功能基因注释

#!/bin/bash
# 脚本名称：13_ncycdb.sh
# 功能：氮循环基因注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - /database/NCycDB/NCycDB.dmnd
# 输出：./3_ncycdb/ncyc_hits.tsv
# 笔记：
# （1）目的：注释氮循环相关功能基因（如 nifH, amoA, nirK 等）。是单步骤功能注释。
# （2）适用性：适用于宏基因组ORFs，包括土壤、水体、肠道等。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：通常用 DIAMOND + NCycDB，可 conda 构建环境。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_ncycdb

diamond blastp \
    --db /database/NCycDB/NCycDB.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_ncycdb/ncyc_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5

echo "✅ NCycDB 注释完成！结果在 ./3_ncycdb/ncyc_hits.tsv"


## 子任务43：使用PCycDB数据库，对磷循环功能基因进行注释

#!/bin/bash
# 脚本名称：17_pcycdb.sh
# 功能：磷循环基因注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - /database/PCycDB/PCycDB.dmnd
# 输出：./3_pcycdb/pcyc_hits.tsv
# 笔记：
# （1）目的：注释磷循环相关基因（如 phoD, phoX, ppk 等）。
# （2）适用性：适用于环境宏基因组（土壤、水体），肠道中较少但可用。
# （3）归类：第三部分
# （4）工具与环境：DIAMOND + PCycDB，可 conda 环境。
# （5）数据准备：proteins.faa

set -e

mkdir -p ./3_pcycdb

diamond blastp \
    --db /database/PCycDB/PCycDB.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_pcycdb/pcyc_hits.tsv \
    --outfmt 6 qseqid sseqid pident evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5

echo "✅ PCycDB 注释完成！结果在 ./3_pcycdb/pcyc_hits.tsv"


## 子任务44：使用MCycDB数据库，分析微生物组的甲烷循环过程

#!/bin/bash
# 脚本名称：18_mcycdb.sh
# 功能：甲烷循环基因注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 数据库：
#   - /database/MCycDB/MCycDB.dmnd
# 输出：./3_mcycdb/mcyc_hits.tsv
# 笔记：
# （1）目的：注释甲烷生成/氧化相关基因（如 mcrA, pmoA 等）。
# （2）适用性：适用于厌氧环境（湿地、瘤胃、沉积物），肠道中较少。
# （3）归类：第三部分
# （4）工具与环境：DIAMOND + MCycDB，conda 可行。
# （5）数据准备：proteins.faa

set -e

mkdir -p ./3_mcycdb

diamond blastp \
    --db /database/MCycDB/MCycDB.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_mcycdb/mcyc_hits.tsv \
    --outfmt 6 qseqid sseqid pident evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5

echo "✅ MCycDB 注释完成！结果在 ./3_mcycdb/mcyc_hits.tsv"


## 子任务48：使用Swiss-Prot数据库，对蛋白质进行注释

#!/bin/bash
# 脚本名称：22_swissprot_annotation.sh
# 功能：基于 DIAMOND 对 proteins.faa 进行 Swiss-Prot 注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# Swiss-Prot DIAMOND 数据库路径：
#   - /database/uniprot/swissprot.dmnd （需提前用 diamond makedb 构建）
# 输出目录：./3_swissprot/
# 输出文件：./3_swissprot/swissprot_hits.tsv
# 笔记：
# （1）目的：利用高可信度、人工审阅的 Swiss-Prot 蛋白数据库对宏基因组预测的蛋白序列进行功能注释，提升注释准确性。这是一个单步骤功能注释任务。
# （2）适用性：适用于所有宏基因组学研究场景，包括微生物组和肠道微生物宏基因组。Swiss-Prot 虽覆盖度低于 NR，但注释质量高，适合高置信度分析。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：主要工具为 diamond（用于快速比对）。可通过 conda 构建独立环境（如 mamba create -n diamond-env; mamba install -c bioconda diamond），在 Linux 中高效完成。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_swissprot

diamond blastp \
    --db /database/uniprot/swissprot.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_swissprot/swissprot_hits.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5 \
    --max-target-seqs 1

echo "✅ Swiss-Prot 注释完成！结果保存至 ./3_swissprot/swissprot_hits.tsv"


## 子任务42：使用Pfam数据库对蛋白质进行结构域层面的注释

#!/bin/bash
# 脚本名称：16_pfam_hmmscan.sh
# 功能：使用 hmmscan 对 proteins.faa 进行 Pfam 结构域注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# Pfam HMM 数据库路径：
#   - /database/Pfam/Pfam-A.hmm （需已用 hmmpress 处理）
# 输出目录：./3_pfam/
# 输出文件：./3_pfam/pfam_domains.tbl
# 笔记：
# （1）目的：识别蛋白质中的保守结构域（domains）和功能位点，从结构层面理解蛋白功能。这是单步骤结构域注释任务。
# （2）适用性：广泛适用于宏基因组蛋白功能解析，尤其在区分同源基因（如 pmoA/amoA）时具有优势，适用于肠道等复杂微生物组。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：推荐使用 hmmscan（来自 HMMER 套件）配合 Pfam HMM 库。可通过 conda 安装（mamba install -c bioconda hmmer pfam），在 Linux 中运行。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_pfam

hmmscan \
    --cpu 32 \
    --tblout ./3_pfam/pfam_domains.tbl \
    --noali \
    /database/Pfam/Pfam-A.hmm \
    ./D_6_ORFs/proteins.faa > /dev/null

echo "✅ Pfam 结构域注释完成！结果保存至 ./3_pfam/pfam_domains.tbl"


## 子任务49：使用NR数据库，对蛋白质进行注释

#!/bin/bash
# 脚本名称：23_nr_annotation.sh
# 功能：基于 DIAMOND 对 proteins.faa 进行 NR 注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# NR DIAMOND 数据库路径：
#   - /database/ncbi/nr.dmnd （需提前构建）
# 输出目录：./3_nr/
# 输出文件：./3_nr/nr_hits.tsv
# 笔记：
# （1）目的：利用 NCBI 非冗余蛋白数据库（NR）对宏基因组蛋白进行最广泛的同源比对与功能推断。这是标准的单步骤注释任务。
# （2）适用性：适用于所有宏基因组研究，因其覆盖度高，常用于初步功能筛查，包括肠道微生物组。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：使用 diamond 比对 NR 数据库。NR 数据库体积大（>200 GB），但可在 Linux 通过 conda 环境运行 DIAMOND。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_nr

diamond blastp \
    --db /database/ncbi/nr.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_nr/nr_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5 \
    --max-target-seqs 1

echo "✅ NR 注释完成！结果保存至 ./3_nr/nr_hits.tsv"


## 子任务41：使用InterProScan对蛋白进行注释和功能分析

#!/bin/bash
# 脚本名称：15_interproscan.sh
# 功能：InterProScan 全面蛋白注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# 输出：./3_interpro/proteins.tsv (TSV格式)
# 笔记：
# （1）目的：通过整合多个数据库（Pfam, Gene3D, SUPERFAMILY等）进行蛋白结构域与功能注释。
# （2）适用性：通用于所有宏基因组蛋白序列。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：interproscan，需 Java，可 conda 安装（但常需手动配置）。
# （5）数据准备：proteins.faa


set -e

mkdir -p ./3_interpro

interproscan.sh \
    -i ./D_6_ORFs/proteins.faa \
    -f TSV \
    -o ./3_interpro/proteins.tsv \
    -cpu 32 \
    -appl Pfam,SMART,PROSITE,PRINTS,Gene3D,SUPERFAMILY

echo "✅ InterProScan 完成！结果在 ./3_interpro/proteins.tsv"


## 子任务51：使用OrthoDB数据库对蛋白质进行注释

#!/bin/bash
# 脚本名称：26_orthodb_annotation.sh
# 功能：基于 DIAMOND 对 proteins.faa 进行 OrthoDB 注释
# 输入：
#   - ./D_6_ORFs/proteins.faa
# OrthoDB DIAMOND 数据库路径：
#   - /database/OrthoDB/odb11v0.dmnd （基于 OrthoDB v11）
# 输出目录：./3_orthodb/
# 输出文件：./3_orthodb/orthodb_hits.tsv
# 笔记：
# （1）目的：将蛋白序列映射到 OrthoDB 的直系同源群（orthologous groups），用于进化与功能保守性分析。这是单步骤注释任务。
# （2）适用性：适用于跨物种比较的宏基因组研究，包括肠道菌群的系统发育功能分析。
# （3）归类：第三部分：基于contigs、ORF层面的分析
# （4）工具与环境：使用 diamond 比对 OrthoDB 蛋白库。OrthoDB 提供预构建的 FASTA，可自行用 diamond makedb 构建索引。conda 环境完全支持。
# （5）数据准备：ORFs（proteins.faa）

set -e

mkdir -p ./3_orthodb

diamond blastp \
    --db /database/OrthoDB/odb11v0.dmnd \
    --query ./D_6_ORFs/proteins.faa \
    --out ./3_orthodb/orthodb_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads 32 \
    --evalue 1e-5 \
    --max-target-seqs 1

echo "✅ OrthoDB 注释完成！结果保存至 ./3_orthodb/orthodb_hits.tsv"


## 子任务16：自定义数据库的比对解释（归类到第三部分其他独立数据库注释子任务之后）
#!/bin/bash
# 脚本名称：9_custom_degradation.sh
# 功能：基于自定义数据库，注释重金属抗性（如 czcA, arsC）和农药降解基因（如 linA, opd）
# 前置数据准备：
#   - proteins.faa
#   - 自定义数据库 custom_degradative.faa（含已知降解酶序列）
#   - 已构建 DIAMOND 数据库：custom_degradative.dmnd
# 输出：
#   - output_dir/degradative_hits.tsv

set -e

output_dir="output_dir"
mkdir -p "${output_dir}"

diamond blastp \
  --db custom_degradative.dmnd \
  --query proteins.faa \
  --out "${output_dir}/degradative_hits.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --threads 16 \
  --evalue 1e-10 \
  --id 70 \
  --query-cover 70

echo "✅ 功能注释完成"


## 子任务52：使用metacompare2.0进行抗性组风险评估
# 由于是进行风险评估，跟物种没有什么关系，因此这一步是在contigs水平进行的。然后，根据样本的覆盖度进行回帖，进而得到单样本的风险，然后再组间求均值。

# 准备工作
# 三部曲配置环境(存在版本冲突,用下面的办法先测试版本)
mamba create -n metacompare2.0 python=3.12
conda activate metacompare2.0
mamba install -c anaconda numpy
mamba install -c anaconda pandas
mamba install -c bioconda biopython
mamba install -c bioconda pprodigal
mamba install -c bioconda diamond=0.9.14
mamba install -c bioconda mmseqs2


# 设置 channel 优先级（测试可以共存的版本,尤其是确定python版本,然后把版本号写在上面）
mamba create -n metacompare2.0 \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  python \
  numpy \
  pandas \
  biopython \
  pprodigal \
  "diamond=0.9.14" \
  mmseqs2 \
  --override-channels \
  --dry-run
  

# 程序+数据库下载（要求：数据库放在与程序同文件夹中）
cd /database/work/zryan/software
git clone https://github.com/mrumi/MetaCompare2.0.git
cd ./MetaCompare2.0
wget https://zenodo.org/api/records/10626079/files/metacmpDB.tar.gz/content
tar -zxvf metacmpDB.tar.gz

# 运行工具

#!/bin/bash
# 脚本名称：C_4_metacompare2.sh
# 功能：使用metacompare2.0进行抗性组风险评估。
# 前置数据准备：
#   - xxx
# 输出：
#   - $OUTPUT_DIR/contig.metacompare.resistomerisk.tsv


# 设置 MetaCompare 2.0 脚本路径
METACOMPARE_SCRIPT="/database/work/zryan/software/MetaCompare2.0/metacompare.py"
# 设置输入 contigs 路径
CONTIGS_FILE="./D_4_contigs/mg_6.contigs.fa"
# 设置输出目录（将自动创建）
OUTPUT_DIR="/database/work/zryan/analysis/mg_6_new/C_4_ResistomeRisk"
mkdir -p "$OUTPUT_DIR"
# 设置线程数
THREADS=72
# 设置风险评分模式：0 = 同时计算生态与人类健康风险
RISK_MODE=0



# 运行 MetaCompare 2.0
python "$METACOMPARE_SCRIPT" \
  -c "$CONTIGS_FILE" \
  -t "$THREADS" \
  -b "$RISK_MODE" \
  -o "$OUTPUT_DIR"

echo "MetaCompare 2.0 analysis completed!"
echo "Output saved to: $OUTPUT_DIR"


# 第四部分：基于MAGs层面的分析
## 子任务17：MAGs丰度定量

#!/bin/bash
# 脚本名称：14_coverm_quant_MAGs.sh
# 功能：使用 CoverM 全面定量弹尾虫 Folsomia candida 肠道微生物 MAGs 在各样本中的多种覆盖度指标
# 前置数据准备：
#   - MAGs 目录: ./D_5_MAGs/（每个 MAG 为一个 .fa 或 .fasta 文件，需要专门用参数指定fa扩展名，否则会按fna识别而出错。）
#   - clean reads: ./D_2_HRreads/cleanreads/ 下的 *.hrm.1.fastq 和 *.hrm.2.fastq（未压缩）
# 输出：
#   - ./M_1_MAGs_quant/coverm_output/ 包含多种覆盖度矩阵

set -e

# 定义路径
CLEAN_READS_DIR="./D_2_HRreads/cleanreads"
MAGS_DIR="./D_5_MAGs"
OUTPUT_DIR="/database_new/work/zryan/TEMP_analysis/mg_6_M_1_MAGs_quant/coverm_output"
TMPDIR="/database_new/work/zryan/TEMP_analysis/mg_6_M_1_MAGs_quant/tmp"

# 创建输出目录与临时目录（确保存在）
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMPDIR}/allmethods"
mkdir -p "${TMPDIR}/coverage_histogram"
mkdir -p "${TMPDIR}/length"


# 运行 coverm genome：length方法（MAGs的--min-covered-fraction默认值为10，length方法要求为0）
# 为了保险，先运行这个，且先设置参数，否则还是会触发报错。
TMPDIR="${TMPDIR}/length" coverm genome \
  --genome-fasta-directory "${MAGS_DIR}" \
  --genome-fasta-extension fa \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --min-covered-fraction 0 \
  --methods length \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT_length.tsv" \


# 运行 coverm genome：coverage_histogram 方法
TMPDIR="${TMPDIR}/coverage_histogram" coverm genome \
  --genome-fasta-directory "${MAGS_DIR}" \
  --genome-fasta-extension fa \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --methods coverage_histogram \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT_coverage_histogram.tsv" \


# 运行 coverm genome：主方法集合（兼容的常规方法）
# 注意：genome 模式下 --methods 不包含 metabat；metabat 是 contig 模式专属

TMPDIR="${TMPDIR}/allmethods" coverm genome \
  --genome-fasta-directory "${MAGS_DIR}" \
  --genome-fasta-extension fa \
  --coupled \
  ${CLEAN_READS_DIR}/CK_1.hrm.1.fastq ${CLEAN_READS_DIR}/CK_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_2.hrm.1.fastq ${CLEAN_READS_DIR}/CK_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_3.hrm.1.fastq ${CLEAN_READS_DIR}/CK_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_4.hrm.1.fastq ${CLEAN_READS_DIR}/CK_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_5.hrm.1.fastq ${CLEAN_READS_DIR}/CK_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/CK_6.hrm.1.fastq ${CLEAN_READS_DIR}/CK_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_1.hrm.1.fastq ${CLEAN_READS_DIR}/T1_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_2.hrm.1.fastq ${CLEAN_READS_DIR}/T1_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_3.hrm.1.fastq ${CLEAN_READS_DIR}/T1_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_4.hrm.1.fastq ${CLEAN_READS_DIR}/T1_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_5.hrm.1.fastq ${CLEAN_READS_DIR}/T1_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T1_6.hrm.1.fastq ${CLEAN_READS_DIR}/T1_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_1.hrm.1.fastq ${CLEAN_READS_DIR}/T3_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_2.hrm.1.fastq ${CLEAN_READS_DIR}/T3_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_3.hrm.1.fastq ${CLEAN_READS_DIR}/T3_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_4.hrm.1.fastq ${CLEAN_READS_DIR}/T3_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_5.hrm.1.fastq ${CLEAN_READS_DIR}/T3_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T3_6.hrm.1.fastq ${CLEAN_READS_DIR}/T3_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_1.hrm.1.fastq ${CLEAN_READS_DIR}/T5_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_2.hrm.1.fastq ${CLEAN_READS_DIR}/T5_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_3.hrm.1.fastq ${CLEAN_READS_DIR}/T5_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_4.hrm.1.fastq ${CLEAN_READS_DIR}/T5_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_5.hrm.1.fastq ${CLEAN_READS_DIR}/T5_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T5_6.hrm.1.fastq ${CLEAN_READS_DIR}/T5_6.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_1.hrm.1.fastq ${CLEAN_READS_DIR}/T6_1.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_2.hrm.1.fastq ${CLEAN_READS_DIR}/T6_2.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_3.hrm.1.fastq ${CLEAN_READS_DIR}/T6_3.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_4.hrm.1.fastq ${CLEAN_READS_DIR}/T6_4.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_5.hrm.1.fastq ${CLEAN_READS_DIR}/T6_5.hrm.2.fastq \
  ${CLEAN_READS_DIR}/T6_6.hrm.1.fastq ${CLEAN_READS_DIR}/T6_6.hrm.2.fastq \
  --methods relative_abundance mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
  --min-read-aligned-length 50 \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 80 \
  --threads 72 \
  --output-file "${OUTPUT_DIR}/STDOUT.tsv" \



# 注：genome 模式不支持 --methods metabat（metabat 是 contig-level 的覆盖度估算方法）

echo "✅ CoverM MAGs 覆盖度定量完成！结果在 ${OUTPUT_DIR}/"
echo "💡 输出包含以下指标：mean, covered_fraction, covered_bases, variance, length, count, reads_per_base, rpkm, tpm, coverage_histogram"

# 说明：参数设置原因（与contigs水平脚本一致，适用于土壤动物肠道宏基因组）
# --min-read-aligned-length 50 \
# 宏基因组 reads 通常为 150 bp（PE），50 bp 是可靠比对的下限。过低会引入噪声；过高会丢失短 contig 的覆盖信息。
# 类似研究常用 50–75 bp（Sczyrba et al., 2017, *Nature Methods*; https://doi.org/10.1038/nmeth.4458）。
#
# --min-read-percent-identity 90 \
# 宏基因组中物种多样性高，90–95% identity 是区分近缘菌株的常用阈值。
# 90% 平衡了灵敏度与特异性，避免将 reads 错配到远缘同源序列（Nayfach et al., 2016, *Genome Research*; https://doi.org/10.1101/gr.201868.115）。
#
# --min-read-aligned-percent 80 \
# 要求至少 80% 的 read 被比对上，防止部分比对（partial alignment）导致的假阳性覆盖。
# 这在短 contigs（<1 kb）中尤为重要（Bishara et al., 2018, *Nature Biotechnology*; https://doi.org/10.1038/nbt.4075）。




## 子任务18/子任务34：使用GTDB-tk工具物种注释与建树（基于 filtered MAGs）（与34相同）
## 子任务18、物种注释与建树（基于MAGs）
# 建议先使用checkM进行质控。参数设置：completeness > 70%, contamination < 10%
# 注意：GTDB-Tk 需要预下载数据库（gtdbtk data），且 de_novo_wf 会构建新树。
# gtdbtk官网：(按照官网上的教程，三部曲+download-db.sh全自动下载对应版本的数据库到gtdbtk的conda环境env文件夹。这样卸载env对应版本，数据库也会删除。)
# https://ecogenomics.github.io/GTDBTk/index.html

# 部署gtdbtk
# 版本需要到官网https://ecogenomics.github.io/GTDBTk/index.html上找最新版本号，不要用2.6.1
mamba create -n gtdbtk-123.456 -c conda-forge -c bioconda gtdbtk=123.456
download-db.sh
conda deactivate
conda activate gtdbtk-123.456 # 重新打开环境才能让环境变量生效。
gtdbtk check_install # 全部OK就可以


# 正式脚本：针对细菌进行分类与建树。
#!/bin/bash
# 脚本名称：15_gtdbtk_classify_and_tree.sh
# 功能：对 MAGs 进行 GTDB 物种注释并构建系统发育树
# 前置数据准备：
#   - filtered_mag_directory/ （MAGs，*.fa）
#   - GTDB-Tk 数据库已配置（通过 gtdbtk data）
# 输出：
#   - classify_output/
#   - tree_output/

set -e

# 1. 物种分类注释
gtdbtk classify_wf \
  --genome_dir filtered_mag_directory \
  --out_dir classify_output \
  --extension fa \
  --cpus 72 \
  --prefix mag

# 2. 从头构建系统发育树（仅细菌,gtdbtk不支持真菌）
gtdbtk de_novo_wf \
  --genome_dir filtered_mag_directory \
  --bacteria \
  --outgroup_taxon p__Patescibacteria \
  --out_dir tree_output \
  --cpus 72 \
  --prefix mag_bacteria

# 3. 构建古菌树
gtdbtk de_novo_wf \
  --genome_dir filtered_mag_directory \
  --archaea \
  --outgroup_taxon p__Altiarchaeota \
  --out_dir tree_output \
  --cpus 72 \
  --prefix mag_archaea

echo "✅ GTDB-Tk 分析完成！分类结果在 classify_output/，树在 tree_output/"

# 扩展版本：针对细菌、古菌进行分类和建树
## 子任务18/34：使用GTDB-tk对MAGs进行物种注释，并分别构建细菌/古菌系统发育树
# 输入：./D_9_filteredMAGs/filtered_mags/*.fa
# 输出：
#   - ./M_2_trees_gtdbtk/classify_output/
#   - ./M_2_trees_gtdbtk/tree_bacteria_output/
#   - ./M_2_trees_gtdbtk/tree_archaea_output/
# 临时文件：./M_2_trees_gtdbtk/tmp/

#!/bin/bash
set -e

# 路径配置
GENOME_DIR="./D_9_filteredMAGs/filtered_mags"
OUTPUT_BASE="./M_2_trees_gtdbtk"
CLASSIFY_OUT="${OUTPUT_BASE}/classify_output"
TREE_BACT_OUT="${OUTPUT_BASE}/tree_bacteria_output"
TREE_ARCH_OUT="${OUTPUT_BASE}/tree_archaea_output"
TMP_DIR="${OUTPUT_BASE}/tmp"
TMP_BACT="${TMP_DIR}/tmp_bacteria"
TMP_ARCH="${TMP_DIR}/tmp_archaea"

EXT="fa"
CPUS=72
PREFIX="mag"

# 创建输出和临时目录
mkdir -p "$CLASSIFY_OUT"
mkdir -p "$TREE_BACT_OUT"
mkdir -p "$TREE_ARCH_OUT"
mkdir -p "$TMP_BACT" "$TMP_ARCH"

# Step 1: 全局分类注释（获取域信息）
echo "🔬 正在运行 GTDB-Tk classify_wf 获取分类信息..."
gtdbtk classify_wf \
  --genome_dir "$GENOME_DIR" \
  --out_dir "$CLASSIFY_OUT" \
  --extension "$EXT" \
  --cpus "$CPUS" \
  --prefix "$PREFIX"

# Step 2: 按域分离MAGs（使用软链接到临时目录）
CLASSIFY_TSV="${CLASSIFY_OUT}/${PREFIX}.gtdbtk.classify_output.tsv"

if [ ! -f "$CLASSIFY_TSV" ]; then
  echo "❌ 错误：分类结果文件不存在: $CLASSIFY_TSV"
  exit 1
fi

# 跳过表头，按域建立软链接
while IFS=$'\t' read -r genome_id classification _; do
  if [[ "$classification" == d__Bacteria* ]]; then
    ln -sf "$(realpath "${GENOME_DIR}/${genome_id}.${EXT}")" "${TMP_BACT}/${genome_id}.${EXT}"
  elif [[ "$classification" == d__Archaea* ]]; then
    ln -sf "$(realpath "${GENOME_DIR}/${genome_id}.${EXT}")" "${TMP_ARCH}/${genome_id}.${EXT}"
  else
    echo "⚠️  警告：未识别域的基因组: $genome_id ($classification)"
  fi
done < <(tail -n +2 "$CLASSIFY_TSV")

# Step 3: 分别构建细菌和古菌的 de novo 树

# 细菌树
if [ -n "$(ls -A "$TMP_BACT" 2>/dev/null)" ]; then
  echo "🌳 正在构建细菌系统发育树..."
  gtdbtk de_novo_wf \
    --genome_dir "$TMP_BACT" \
    --bacteria \
    --outgroup_taxon p__Patescibacteriota \
    --out_dir "$TREE_BACT_OUT" \
    --extension "$EXT" \
    --cpus "$CPUS" \
    --prefix "${PREFIX}_bacteria"
else
  echo "ℹ️  无细菌基因组，跳过细菌树构建。"
fi

# 古菌树
if [ -n "$(ls -A "$TMP_ARCH" 2>/dev/null)" ]; then
  echo "🌳 正在构建古菌系统发育树..."
  gtdbtk de_novo_wf \
    --genome_dir "$TMP_ARCH" \
    --archaea \
    --outgroup_taxon p__Altiarchaeota \
    --out_dir "$TREE_ARCH_OUT" \
    --extension "$EXT" \
    --cpus "$CPUS" \
    --prefix "${PREFIX}_archaea"
else
  echo "ℹ️  无古菌基因组，跳过古菌树构建。"
fi

# Step 4: 清理临时目录（可选：如需保留中间文件用于调试，可注释掉此行）
rm -rf "$TMP_DIR"

echo "✅ GTDB-Tk 分析完成！"
echo "   - 分类结果: $CLASSIFY_OUT/"
echo "   - 细菌树:   $TREE_BACT_OUT/"
echo "   - 古菌树:   $TREE_ARCH_OUT/"



## 子任务19：可移动遗传原件MGEs功能注释（基于MAGs）
# 📌 流程位置与前置要求
# 位置：在获得 MAGs 的蛋白质序列 之后（通常由 Prodigal 预测）。
# 前置步骤：
# 安装 mobileOG-db 及其脚本（mobileOGs-pl-kyanite.sh）。
# 不需要将MAGs转化为蛋白质，程序中有。
# 注意：你提供的命令是针对单个 MAG 的，需扩展为多样本。
# 教程:(同样适用于contigs级别)
# https://zhuanlan.zhihu.com/p/688079390
# https://github.com/clb21565/mobileOG-db

# 准备工作
# 配置conda环境
mamba create -n mobileOG-db python=3.6.15
mamba activate mobileOG-db
mamba install -conda-forge biopython
mamba install -c bioconda prodigal
mamba install -c bioconda diamond
mamba install -c anaconda pandas

# 下载mobileOG-db程序
git clone https://github.com/clb21565/mobileOG-db?tab=readme-ov-file
chmod +x ./mobileOG-pl/mobileOGs-pl-kyanite.sh

# 下载数据库（直接去这个网站手动下载，然后ftp上传服务器)
https://mobileogdb.flsi.cloud.vt.edu/entries/database_download

# 把下载的faa蛋白数据库，转化为diamond数据库格式(建库)
cd ./path_to_database 
diamond makedb --in ./mobileOG-db_beatrix-1.6.All.faa -d ./mobileOG-db_beatrix-1.6.All.dmnd

# MGEs标注
#!/bin/bash
# 脚本名称：16_mobileog_mge_annotation.sh
# 功能：使用 mobileOG-db 注释 filtered MAGs 中的可移动遗传元件（MGEs）
# 流程：Prodigal → DIAMOND (main) → Python summary + DIAMOND (gene list)
# 前置要求：
#   - Conda 环境 'mobileOG-db' 已激活（含 prodigal, diamond, python=3.6.15, pandas）
#   - 输入为高质量 MAGs 的核苷酸 FASTA 文件（.fa）
# 输出：
#   - 每个 MAG 对应子目录：含 .tsv, .csv, gene_mobileOG.list 等

set -e

# ==============================
# 📌 路径与参数配置
# ==============================

# 输入 MAGs 目录
INPUT_MAGS_DIR="./D_9_filteredMAGs/filtered_mags"
INPUT_MAGS_PATTERN="$INPUT_MAGS_DIR/*.fa"

# DIAMOND 数据库与元数据
DIAMOND_DB="/database/work/zryan/biodb/mobileOG-db_1.6/mobileOG-db_beatrix-1.6.All.dmnd"
METADATA_CSV="/database/work/zryan/biodb/mobileOG-db_1.6/mobileOG-db-beatrix-1.6-All.csv"

# Python 脚本路径（用于结果汇总）
MOBILEOG_PY="/database/work/zryan/software/mobileOG-db/mobileOG-pl/mobileOGs-pl-kyanite.py"

# DIAMOND 参数（主比对）
KVALUE=15
ESCORE="1e-5"
PIDENTVALUE=90
QUERYSCORE=90
DIAMOND_THREAD=80

# 输出主目录
OUTPUT_DIR="./M_3_MGEs_mobileOGdb"

# ==============================
# 🛠️ 脚本主体
# ==============================

# 检查依赖文件
if [ ! -f "$DIAMOND_DB" ]; then
  echo "❌ 错误：DIAMOND 数据库未找到！$DIAMOND_DB"
  exit 1
fi
if [ ! -f "$METADATA_CSV" ]; then
  echo "❌ 错误：元数据 CSV 未找到！$METADATA_CSV"
  exit 1
fi
if [ ! -f "$MOBILEOG_PY" ]; then
  echo "❌ 错误：Python 脚本未找到！$MOBILEOG_PY"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# 定义单个 MAG 处理函数
annotate_mge() {
  local mag_fasta="$1"
  local mag_fasta_abs=$(readlink -f "$mag_fasta")
  local base_name=$(basename "$mag_fasta" .fa)
  local out_dir="$OUTPUT_DIR/${base_name}_mobileog"

  echo "[$(date)] START: $base_name"
  mkdir -p "$out_dir"

  (
    cd "$out_dir" || exit 1

    # Step 1: Prodigal
    echo "→ Prodigal..."
    prodigal -i "$mag_fasta_abs" -p meta -a "${base_name}.faa" -o /dev/null

    # Step 2: DIAMOND (main)
    echo "→ DIAMOND (main for summary)..."
    diamond blastp \
      -q "${base_name}.faa" \
      --db "$DIAMOND_DB" \
      --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend \
      -k "$KVALUE" \
      -o "${base_name}.tsv" \
      -e "$ESCORE" \
      --query-cover "$QUERYSCORE" \
      --id "$PIDENTVALUE" \
      --threads "$DIAMOND_THREAD" \
      --ultra-sensitive

    # Step 3: Python summary —— 关键：--o 只传 base_name！
    echo "→ Python summary..."
    python "$MOBILEOG_PY" \
      --i "${base_name}.tsv" \
      --o "$base_name" \
      -m "$METADATA_CSV"

    # Step 4: DIAMOND for gene list
    echo "→ DIAMOND (for gene list)..."
    diamond blastp \
      --db "$DIAMOND_DB" \
      --query "${base_name}.faa" \
      --outfmt 6 \
      --threads "$DIAMOND_THREAD" \
      --max-target-seqs 1 \
      -e "$ESCORE" \
      --ultra-sensitive \
      --out "${base_name}_gene_diamond.f6"

    # Step 5: gene list
    echo "→ Generating gene_mobileOG.list..."
    cut -f 1,2 "${base_name}_gene_diamond.f6" | uniq | \
      sed '1 i Name\tResGeneID' > "${base_name}_gene_mobileOG.list"

    echo "→ DONE: $base_name"
  )

  echo "[$(date)] FINISH: $base_name"
}

# 获取 MAG 列表
mag_files=$(ls $INPUT_MAGS_PATTERN 2>/dev/null)
if [ -z "$mag_files" ]; then
  echo "❌ 错误：未在 $INPUT_MAGS_DIR 中找到 .fa 文件！"
  exit 1
fi

# 单线程处理
echo "🔍 找到 $(echo $mag_files | wc -w) 个 MAGs，开始单线程注释..."
for mag in $mag_files; do
  echo "→ Processing: $(basename "$mag")"
  annotate_mge "$mag"
done

echo "✅ mobileOG MGE 注释完成！结果保存在 $OUTPUT_DIR/"


## 子任务20：使用 DIAMOND 映射 MAGs 蛋白质到 ARGs 数据库（如 SARG）
# 📌 流程位置与前置要求
# 位置:MAGs 的蛋白质序列（.faa）。
# 前置步骤：
# Prodigal 预测 MAG 蛋白质 → mag_proteins/*.faa；
# 已构建 SARG 或 CARD 等 ARG 数据库的 DIAMOND 索引（.dmnd）。
# 注意：你提供的命令是单样本，需改为多样本。

#!/bin/bash
# 脚本名称：17_diamond_arg_annotation.sh
# 功能：使用 DIAMOND 将 MAG 蛋白质比对到 ARG 数据库（如 SARG）
# 前置数据准备：
#   - mag_proteins/ （MAG 蛋白质，*.faa）
#   - sarg_database.dmnd （已构建的 DIAMOND 数据库）
# 输出：
#   - arg_results/ （每个 MAG 一个 .txt 结果文件）

set -e

mkdir -p arg_results

# 定义单样本注释函数
blast_arg() {
  local prot_file=$1
  local base_name=$(basename "$prot_file" .faa)
  local out_file="arg_results/${base_name}_arg.txt"

  diamond blastp \
    --threads 12 \
    -d sarg_database.dmnd \
    -q "$prot_file" \
    -e 1e-10 \
    --id 90 \
    --query-cover 90 \
    --max-target-seqs 1 \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    -o "$out_file"
}

export -f blast_arg

prot_files=$(ls mag_proteins/*.faa)

# 并行运行（-j 10，每任务12线程，总约120线程，符合你服务器72核+超线程）
parallel -j 10 blast_arg ::: $prot_files

echo "✅ ARG 注释完成！结果在 arg_results/"


## 子任务24/子任务33：使用CarveMe进行微生物群落基因组尺度代谢建模（与33相同）
# 工具：CarveMe、ModelSEED
# 输入：MAGs 的 KO 注释
# 输出：基因组尺度代谢模型（GEMs），预测碳源利用、代谢互作

# 这是2024年6月发表在《Trends in Endocrinology & Metabolism》上发表的名为“Emerging methods for genome-scale metabolic modeling of microbial communities”的研究。
# 揭秘微生物社区的“社交网络”：基因组规模代谢模型(GEMs)的新前沿


### 数据准备
# 目的：得到KO注释。

#!/bin/bash
# 脚本名称：23_prepare_ko_annotations.sh
# 功能：从 KEGG 注释结果（如 8_kegg_annotation.sh 的 ko_output.txt）生成 per-MAG KO 列表
# 前置数据准备：
#   - mag_proteins/（MAG 蛋白质 .faa）
#   - ko_output.txt（来自 exec_annotation，格式：protein_id<tab>Kxxxxx）
# 输出：
#   - ko_annotations/（每个 MAG 一个 .ko 文件）

mkdir -p ko_annotations

# 假设 ko_output.txt 格式为：prot_id<tab>KO_id
while read prot_file; do
  base_name=$(basename "$prot_file" .faa)
  # 提取该 MAG 所有蛋白 ID 前缀（假设蛋白 ID 为 bin.1_00001）
  # 此处简化：假设蛋白 ID 以 ${base_name}_ 开头
  awk -v prefix="${base_name}_" '$1 ~ "^"prefix {print $1 "\t" $2}' ko_output.txt > "ko_annotations/${base_name}.ko"
done < <(ls mag_proteins/*.faa)


### GEMs构建

#!/bin/bash
# 脚本名称：23_gem_modeling.sh
# 功能：基于 MAGs 的 KO 注释构建基因组尺度代谢模型（GEMs）
# 前置数据准备：
#   - filtered_mags/（高质量 MAGs，*.fa）
#   - ko_annotations/（每个 MAG 对应一个 .ko 文件，格式：gene_id<tab>Kxxxxx）
#     → 若无，需先运行 KEGG 注释（如 8_kegg_annotation.sh）并转换为 KO 列表
# 输出：
#   - output_dir/carveme_models/（.xml 格式 GEMs）

set -e

output_dir="output_dir"
mkdir -p "${output_dir}/carveme_models"

# 检查 KO 注释是否存在
if [ ! -d "ko_annotations" ]; then
  echo "❌ 错误：ko_annotations/ 目录不存在！请先运行 KEGG 注释并将结果转换为每 MAG 一个 .ko 文件。"
  exit 1
fi

# 定义单 MAG 建模函数
build_gem() {
  local mag_file=$1
  local base_name=$(basename "$mag_file" .fa)
  local ko_file="ko_annotations/${base_name}.ko"
  local model_out="${output_dir}/carveme_models/${base_name}.xml"

  if [ ! -f "$ko_file" ]; then
    echo "⚠️  Warning: KO file missing for ${base_name}, skipping."
    return
  fi

  # 使用 CarveMe 构建模型（默认细菌模式）
  carve --fbc2 --gapseq --output "${model_out}" "$ko_file"
}

export -f build_gem

# 并行构建
mag_files=$(ls filtered_mags/*.fa)
parallel -j 8 build_gem ::: $mag_files

echo "✅ 代谢模型构建完成！GEMs 在 ${output_dir}/carveme_models/"



## 子任务25：群体基因组学（within-species variation）
# 如果要进行这个分析,还需要去补一下知识.因为这部分不一定适用于微生物群落。看到的文献都是单物种的群体基因组学。
# 方法：将 reads 回贴到高质量 MAGs，用 Snippy 或 iVar call SNP
# 输出：SNP 密度、dN/dS、选择压力分析


#!/bin/bash
# 脚本名称：24_within_species_snp.sh
# 功能：将 reads 回贴到高质量 MAGs，call SNP 并计算 dN/dS
# 前置数据准备：
#   - filtered_mags/（高质量 MAGs，*.fa）
#   - ./clean_reads/ 下的 clean reads (.R1.qualified.fastq.gz 等)
#   - Snippy 已安装
# 输出：
#   - output_dir/snippy_results/（每个样本-MAG 对的 SNP 结果）
#   - output_dir/dnds_summary.tsv

set -e

output_dir="output_dir"
tmp_dir="tmp_dir"
mkdir -p "${output_dir}/snippy_results" "${tmp_dir}"

# 构建 MAGs 合并索引（Snippy 支持多参考，但通常 per-MAG 更清晰）
# 此处按每个 MAG 单独处理

samples=$(ls ./clean_reads/*.R1.qualified.fastq.gz | xargs -n1 basename | sed 's/\.R1\.qualified\.fastq\.gz//')
mag_files=$(ls filtered_mags/*.fa)

# 定义 SNP calling 函数
call_snp_per_mag_sample() {
  local mag_file=$1
  local sample=$2
  local out_dir="${output_dir}/snippy_results/${sample}_$(basename "$mag_file" .fa)"

  snippy \
    --cpus 8 \
    --outdir "$out_dir" \
    --ref "$mag_file" \
    --R1 "./clean_reads/${sample}.R1.qualified.fastq.gz" \
    --R2 "./clean_reads/${sample}.R2.qualified.fastq.gz"
}

export -f call_snp_per_mag_sample

# 生成所有 (sample, mag) 组合
for mag in $mag_files; do
  for s in $samples; do
    echo "$mag $s"
  done
done | parallel -j 6 call_snp_per_mag_sample {1} {2}

# 后续可使用 snp-sites 或 custom script 计算 dN/dS（此处略）

echo "✅ SNP calling 完成！结果在 ${output_dir}/snippy_results/"





# 第五部分：同时进行多层次的分析
## 子任务21：水平基因转移（HGT）检测
# 工具：DeepVirFinder（病毒）、MobileElementFinder、ISEScan（插入序列）、ICEfinder（整合性接合元件）
# 输入：MAGs 或 contigs
# 意义：解析 ARGs/VFs 的传播机制（如 plasmid vs chromosome）
# 多用途说明：该脚本既可以在contigs层面运行，也可以在MAGs层面运行。若需在 MAG 水平运行（如区分染色体 vs 质粒 ARGs），建议先运行 18_prodigal_mag_to_proteins.sh 和 19_checkm_filter_mags.sh，然后对 filtered_mags/ 中的每个 .fa 文件单独运行上述工具。但为保持流程简洁，本脚本默认在 contigs 水平运行（覆盖更广）。

#!/bin/bash
# 脚本名称：20_hgt_detection.sh
# 功能：对 contigs 或 MAGs 进行 HGT 相关元件检测（IS、ICE、MGE、病毒）
# 前置数据准备：
#   - contigs.fa（用于全组装水平 HGT 检测）
#   - mag_directory/（可选，若需在 MAG 水平运行）
#   - ISEScan、ICEfinder、MobileElementFinder、DeepVirFinder 已安装并配置好数据库
# 输出：
#   - output_dir/isescan/
#   - output_dir/icefinder/
#   - output_dir/mobileelementfinder/
#   - output_dir/deepvirfinder/

set -e

output_dir="output_dir"
tmp_dir="tmp_dir"
mkdir -p "${output_dir}/isescan" "${output_dir}/icefinder" "${output_dir}/mobileelementfinder" "${output_dir}/deepvirfinder" "${tmp_dir}"

# === 1. ISEScan: 插入序列检测 ===
echo "🔍 Running ISEScan on contigs.fa..."
isescan.py contigs.fa /path/to/isescan_db "${output_dir}/isescan" --threads 16

# === 2. ICEfinder: 整合性接合元件检测 ===
echo "🧬 Running ICEfinder..."
ICEfinder.py -i contigs.fa -o "${output_dir}/icefinder" -t 16

# === 3. MobileElementFinder: 可移动元件（含质粒、转座子等）===
echo "🚚 Running MobileElementFinder..."
MobileElementFinder.py predict -i contigs.fa -o "${output_dir}/mobileelementfinder" --db_dir /path/to/mef_db --threads 16

# === 4. DeepVirFinder: 病毒序列识别（用于辅助判断病毒介导 HGT）===
echo "🦠 Running DeepVirFinder..."
deepvirfinder.py -i contigs.fa -o "${output_dir}/deepvirfinder" -l 5000 -m 16

echo "✅ HGT 元件检测完成！结果在 ${output_dir}"

## 子任务35：使用DRAM工具进行微生物组的代谢基因注释（重新写脚本，使得其在一个脚本中，既进行ORFs层面的注释分析，又进行MAGs层面的注释分析）


#!/bin/bash
# 脚本名称：9_dram_annotate.sh
# 功能：DRAM 代谢注释
# 输入：
#   - ./4_filtered_MAGs/fbin.*.fa
# 输出目录：./4_dram_out/
# 笔记：
# （1）目的：对MAGs或基因集进行深度代谢功能注释（含CAZy、KEGG、VFDB等）。是完整注释流程。
# （2）适用性：适用于MAGs或ORFs，广泛用于肠道等微生物组。
# （3）归类：第五部分
# （4）工具与环境：DRAM，conda安装，Linux支持。
# （5）数据准备：filtered MAGs（fbin.*.fa）

set -e

DRAM.py annotate \
    -i ./4_filtered_MAGs \
    -o ./4_dram_out \
    --threads 32

echo "✅ DRAM 注释完成！结果在 ./4_dram_out/"



# 第六部分：独立完整的生物信息学流程
## 子任务32:使用metaGEM工具，从宏基因组层面重建基因组尺度的代谢模型

#!/bin/bash
# 脚本名称：6_metagem_run.sh
# 功能：启动 metaGEM 全流程（从 reads 到 GEMs）
# 输入：
#   - raw reads: ./raw_reads/sample*.R*.raw.fastq.gz
# 配置文件：./metaGEM/config/config.yaml（需用户预先配置样本路径、线程等）
# 输出：MAGs, GEMs, FBA results 等
# 依赖：Snakemake, Conda, metaGEM 已克隆
# 笔记:
# （1）目的：端到端从宏基因组数据重建MAGs并构建GEMs，模拟群落代谢互作。是完整整合流程。
# （2）适用性：专为宏基因组设计，适用于复杂微生物组（含肠道）。
# （3）归类：第六部分：独立完整的生物信息学流程
# （4）工具与环境：metaGEM（Snakemake流程），需 conda 环境，Linux支持。
# （5）数据准备：raw reads（或 clean reads）
# ⚠️ metaGEM 是 Snakemake 流程，需配置 config.yaml。此处提供启动脚本。
# 这里的内容可能有错误，需要到github里面重新弄。

set -e

cd metaGEM

snakemake --use-conda --cores 64 all

echo "✅ metaGEM 流程完成！结果在 ./metaGEM/results/"


## 子任务37：使用MetaWRAP工具，从宏基因组数据挖掘单菌基因组bins（从质控到分箱+注释的独立整体流程）
## 子任务37：使用MetaWRAP工具，从宏基因组数据挖掘单菌基因组bins（从质控到分箱+注释的独立整体流程）

#!/bin/bash
# 脚本名称：11_metawrap_pipeline.sh
# 功能：MetaWRAP 全流程（reads → MAGs）
# 输入：
#   - ./raw_reads/sample*.R*.raw.fastq.gz
# 输出目录：./6_metawrap_out/
# 笔记：
# （1）目的：端到端从 raw reads 到高质量 MAGs（含质控、组装、分箱、去冗余、注释）。是完整独立流程。
# （2）适用性：广泛用于宏基因组MAGs重建，适用于肠道等复杂样本。
# （3）归类：第六部分：独立完整的生物信息学流程
# （4）工具与环境：metawrap，conda 安装，Linux 支持良好。
# （5）数据准备：raw reads（sample*.R1.raw.fastq.gz, sample*.R2.raw.fastq.gz）

set -e

# Step 1: 质控
metawrap read_qc -1 ./raw_reads/*R1.raw.fastq.gz -2 ./raw_reads/*R2.raw.fastq.gz -t 72 -o ./6_metawrap_out/1_qc

# Step 2: 组装
metawrap assembly -1 ./6_metawrap_out/1_qc/*_1.fastq -2 ./6_metawrap_out/1_qc/*_2.fastq -m 800 -t 72 -o ./6_metawrap_out/2_assembly

# Step 3: 分箱（多工具）
metawrap binning -a ./6_metawrap_out/2_assembly/final_assembly.fasta \
    -o ./6_metawrap_out/3_binning \
    --metabat2 --maxbin2 --concoct \
    -1 ./6_metawrap_out/1_qc/*_1.fastq \
    -2 ./6_metawrap_out/1_qc/*_2.fastq \
    -t 72 -m 800

# Step 4: Bin refinement
metawrap bin_refinement -o ./6_metawrap_out/4_refinement \
    -A ./6_metawrap_out/3_binning/metabat2_bins \
    -B ./6_metawrap_out/3_binning/maxbin2_bins \
    -C ./6_metawrap_out/3_binning/concoct_bins \
    -t 32 -m 400

# Step 5: 注释（可选）
metawrap annotate_bins -b ./6_metawrap_out/4_refinement/metawrap_50_10_bins -o ./6_metawrap_out/5_annotation -t 32

echo "✅ MetaWRAP 全流程完成！MAGs 在 ./6_metawrap_out/4_refinement/"



# 第七部分：病毒组专题分析
## 子任务22：病毒组挖掘（Virome）
# 工具：VirSorter2、DeepVirFinder、VIBRANT
# 输入：contigs（>5 kb）或 MAGs
# 输出：病毒 contigs、完整性评分、宿主预测

#!/bin/bash
# 脚本名称：21_viral_contig_mining.sh
# 功能：从 contigs 中挖掘病毒序列，评估完整性并预测宿主
# 前置数据准备：
#   - contigs.fa（建议长度 >1 kb，但 VirSorter2/VIBRANT 内部会过滤）
#   - VirSorter2、VIBRANT、DeepVirFinder 已安装并配置数据库
#   - （可选）高质量 MAGs 目录 filtered_mags/（用于宿主预测）
# 输出：
#   - output_dir/virsorter2/
#   - output_dir/vibrant/
#   - output_dir/viral_summary.tsv（整合结果）

set -e

output_dir="output_dir"
tmp_dir="tmp_dir"
mkdir -p "${output_dir}/virsorter2" "${output_dir}/vibrant" "${tmp_dir}"

# 提取 >5kb contigs（符合多数病毒工具要求）
seqkit seq -m 5000 contigs.fa > "${tmp_dir}/contigs_5k.fa"

# === 1. VirSorter2 ===
echo "🦠 Running VirSorter2..."
virsorter run -w "${output_dir}/virsorter2" -i "${tmp_dir}/contigs_5k.fa" --include-groups dsDNAphage,ssDNA --min-length 5000 -j 16 all

# === 2. VIBRANT ===
echo "🌀 Running VIBRANT..."
VIBRANT_run.py -i "${tmp_dir}/contigs_5k.fa" -folder "${output_dir}/vibrant" -t 16 -virome

# === 3. （可选）整合 DeepVirFinder 结果（若 20_hgt_detection.sh 已运行）
# 此处假设 DeepVirFinder 结果在 ../20_hgt_detection/output_dir/deepvirfinder/
# 用户可根据实际路径调整或重复运行：
# deepvirfinder.py -i "${tmp_dir}/contigs_5k.fa" -o "${output_dir}/deepvirfinder" -l 5000 -m 16

# === 4. 简易整合（示例：提取 VirSorter2 高可信病毒）===
if [ -f "${output_dir}/virsorter2/final-viral-score.tsv" ]; then
  awk -F'\t' '$2 >= 0.9 || ($3 ~ /hallmark/ && $2 >= 0.7)' "${output_dir}/virsorter2/final-viral-score.tsv" | cut -f1 > "${output_dir}/high_confidence_viral_ids.txt"
  seqkit grep -f "${output_dir}/high_confidence_viral_ids.txt" contigs.fa > "${output_dir}/high_confidence_viral_contigs.fa"
fi

echo "✅ 病毒挖掘完成！高置信病毒序列在 ${output_dir}/high_confidence_viral_contigs.fa"


## 子任务23：CRISPER-Cas系统识别
## 子任务23、CRISPER-Cas系统识别
# 工具：PILER-CR、CRISPRCasFinder
# 意义：反映宿主-病毒互作历史，可用于追踪噬菌体感染

#!/bin/bash
# 脚本名称：22_crispr_cas_identification.sh
# 功能：识别 contigs 或 MAGs 中的 CRISPR 阵列与 Cas 基因
# 前置数据准备：
#   - contigs.fa（或 filtered_mags/*.fa）
#   - PILER-CR 和 CRISPRCasFinder 已安装
# 输出：
#   - output_dir/pilercr/
#   - output_dir/crisprcasfinder/

set -e

output_dir="output_dir"
mkdir -p "${output_dir}/pilercr" "${output_dir}/crisprcasfinder"

# === 1. PILER-CR: 快速识别 CRISPR 重复-间隔阵列 ===
echo "🧬 Running PILER-CR..."
pilercr -in contigs.fa -out "${output_dir}/pilercr/crispr_pilercr.txt"

# === 2. CRISPRCasFinder: 识别 CRISPR + Cas 系统（更全面）===
echo "🛡️ Running CRISPRCasFinder..."
CRISPRCasFinder.pl -in contigs.fa -out "${output_dir}/crisprcasfinder" -def General -cas -att -keep

echo "✅ CRISPR-Cas 系统识别完成！结果在 ${output_dir}"


## 子任务36：使用DRAM-v工具进行病毒组代谢基因注释

#!/bin/bash
# 脚本名称：10_dramv_annotate.sh
# 功能：DRAM-v 病毒代谢注释
# 输入：
#   - ./7_viral/viral_bins/*.fa
# 输出目录：./7_dramv_out/
# （1）目的：对病毒基因组进行代谢潜力注释（如辅助代谢基因AMGs）。
# （2）适用性：仅适用于病毒组（vMAGs或病毒contigs）。
# （3）归类：第七部分：病毒组专题分析
# （4）工具与环境：DRAM-v（DRAM的扩展），conda安装。
# （5）数据准备：病毒基因组（viral_contigs.fa 或 vMAGs）

set -e

DRAM-v.py annotate \
    -i ./7_viral/viral_bins \
    -o ./7_dramv_out \
    --threads 32

echo "✅ DRAM-v 注释完成！结果在 ./7_dramv_out/"



# 第八部分：多组学联合分析
## 子任务26：宏基因组与宏转录组整合分析



# 第九部分：扩增子联合分析
## 子任务28：使用Barrnap工具，从宏基因组中提取出16sRNA

#!/bin/bash
# 脚本名称：2_barrnap_16s.sh
# 功能：从 MAGs 中批量提取 16S rRNA 基因
# 输入：
#   - ./4_MAGs/bin.*.fa
# 输出目录：./1_barrnap_16s/
# 输出文件：每个 MAG 对应一个 .16s.fa 文件
# 笔记：
# （1）目的：从组装后的基因组或MAGs中快速识别并提取16S rRNA基因序列。是单步骤注释任务。
# （2）适用性：适用于宏基因组学中的MAGs或参考基因组分析，但不适用于raw reads或contigs未拼接情况；可用于肠道微生物MAGs的16S提取。
# （3）归类：第一部分：数据准备模块（按题设要求）
# （4）工具与环境：barrnap，可通过 conda 安装，在 Linux 中独立运行。
# （5）数据准备：MAGs（bin.1.fa, bin.2.fa...）或 filtered MAGs（fbin.*.fa）



set -e

mkdir -p ./1_barrnap_16s

for mag in ./4_MAGs/bin.*.fa; do
    prefix=$(basename "$mag" .fa)
    barrnap "$mag" --outseq "./1_barrnap_16s/${prefix}.16s.fa" --threads 4
done

echo "✅ Barrnap 完成！16S 序列保存至 ./1_barrnap_16s/"



## 子任务38：SortMeRNA（V4）工具：宏基因组提取rRNA序列

#!/bin/bash
# 脚本名称：12_sortmerna_rRNA.sh
# 功能：从 clean reads 中提取 rRNA 序列
# 输入：
#   - ./1_kneaddata/sample1_hrm_1.fastq, sample1_hrm_2.fastq ...
# 输出目录：./1_sortmerna/
# 输出：sample1_hrm_rRNA.fq（rRNA reads），sample1_hrm_other.fq（非rRNA）
# 笔记：
# （1）目的：从原始reads或contigs中分离rRNA序列（用于后续16S分析或去除rRNA）。是单步骤过滤。
# （2）适用性：适用于宏基因组和宏转录组，常用于数据预处理。
# （3）归类：第九部分：扩增子联合分析
# （4）工具与环境：sortmerna，conda安装，Linux支持。
# （5）数据准备：raw reads 或 clean reads（此处以 clean reads 为例）

set -e

mkdir -p ./1_sortmerna

for r1 in ./1_kneaddata/*_hrm_1.fastq; do
    r2=${r1/_1.fastq/_2.fastq}
    sample=$(basename "$r1" | sed 's/_hrm_1.fastq//')
    
    sortmerna \
        --ref /opt/anaconda3/envs/sortmerna/share/sortmerna/rRNA_databases/*.fasta \
        --reads "$r1" \
        --reads "$r2" \
        --aligned ./1_sortmerna/"${sample}"_hrm_rRNA \
        --other ./1_sortmerna/"${sample}"_hrm_other \
        --fastx \
        --threads 8
done

echo "✅ SortMeRNA 完成！非rRNA reads 保存至 ./1_sortmerna/*_other.fq"

