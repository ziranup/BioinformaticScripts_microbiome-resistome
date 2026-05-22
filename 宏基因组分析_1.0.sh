##########
# 宏基因组分析流程_1.0 #
##########



###### 确定最大运行参数 ######
# 查看逻辑CPU数目
grep -c ^processor /proc/cpuinfo
maxcore
# 我已经用alias把maxcore写进.bashrc了，很方便。
# 输出结果中关键参数：
# CPU(s): 表示逻辑 CPU 核心数（包含超线程，是服务器可同时处理的任务数上限，最关键指标）。
# 例如，输出若显示 CPU(s): 64，表示服务器最多可同时处理 64 个线程任务。
# 计算最大值
# 单个fastp进程用THREADS个线程，同时运行PARALLEL_JOBS个样本，则：
# 总线程数 = THREADS × PARALLEL_JOBS。这个总线程数应 ≤ 服务器逻辑核心数（避免资源过载）。
# 具体计算：
# 假设通过上述命令查到服务器逻辑核心数为 N（例如N=80），则：
# THREADS × PARALLEL_JOBS ≤ N
# 推荐设置逻辑：
# THREADS（单个样本的线程数）：fastp是多线程优化较好的工具，单个进程建议分配 8-16 线程（线程数过高可能因内部调度开销导致效率下降，实测 16 线程是较优值）。最大值建议不超过 24（超过后边际效益递减）。
# PARALLEL_JOBS（同时运行的样本数）：
# 由总核心数和THREADS反推：
# PARALLEL_JOBS_max = 逻辑核心数 ÷ THREADS
# 例如：若逻辑核心数 = 80，THREADS=16，则PARALLEL_JOBS_max=64÷16=5
# 第三步：实际设置建议
# 留有余地：
# 服务器可能同时运行其他任务，建议总线程数设为逻辑核心数的 80%-90%（避免资源耗尽）。
# 例如逻辑核心数 = 64，总线程数可设为 56（而非 64），即16×3=48或14×4=56。
# 测试优化：
# 可先以较小参数（如THREADS=16，PARALLEL_JOBS=2）运行，用top命令观察 CPU 使用率：
#     若 CPU 使用率长期低于 80%，可适当提高PARALLEL_JOBS；
#     若 CPU 使用率接近 100% 且任务卡顿，需降低PARALLEL_JOBS。
# 超线程的影响：
#     若服务器开启超线程（逻辑核心数 = 2× 物理核心数），对于fastp这类计算密集型任务，基于物理核心数计算可能更高效（例如物理核心 32，逻辑核心 64，总线程数可设为 32-48）。




######数据预处理######

# 解压gz文件
gunzip *.gz

# 使用fastp进行数据质控
# 单端质控
fastp -i input_single.fq -o output_single.fq -w 16
# 双端质控（见1_QC.sh)
fastp -i input_R1.fastq -I input_R2.fastq -o output_R1.fastq -O output_R2.fastq -w 16

# 使用bowtie2去宿主
# 建立索引
bowtie2-build host_genome.fa host_index
# 比对及去宿主
# 双端序列比对
bowtie2 --very-sensitive -t -p 24 -x host_index -1 cleaned_R1.fastq -2 cleaned_R2.fastq -S output.sam --un-conc host_removed
# 单端序列比对
bowtie2 --very-sensitive -t -p 24 -x host_index -U cleaned_reads.fq -S output.sam

##备注：常用的基因组地址
# 跳虫内生沃尔巴克氏属:Wolbachia endosymbiont of Folsomia candida
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=169402&mode=Info
# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/169402/
# 跳虫:Folsomia candida
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=158441&mode=Info
# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/158441/




######基于去宿主reads的分析######
######微生物群落与多样性######
# 联用Kracken2和Bracken进行微生物群落分析
# 先用Kraken2,进行物种分类注释。
# 对双端测序数据（sample_R1.fastq/sample_R2.fastq）进行物种水平的分类注释，输出分类结果报告。
# 最终汇总所有读段的分类结果，生成报告（sample.kreport），包含各分类单元（如门、属、种）的读段计数、相对丰度等信息。
kraken2 --threads 40 --db kraken2_database --paired sample_R1.fastq sample_R2.fastq --report sample.kreport

# 然后运行Bracken，进行物种丰度的校正。
# 校正 Kraken2 的分类结果，提高物种水平（Species）的丰度准确性。Kraken2 的分类结果可能存在偏差：短读段可能无法精确匹配到 “物种” 水平（可能匹配到更高分类等级，如属），导致物种丰度被低估。
# 最终输出校正后的物种丰度结果，更接近真实群落组成。
bracken -d kraken2_database -i sample.kreport -o sample.S.bracken -w sample.S.bracken.kreport -l S -t 32

# 使用krakentools合并结果
# krakentools 是一套辅助处理 Kraken/Bracken 结果的脚本集。
# 给予脚本最高权限
# 给combine_kreports.py脚本添加可执行权限（后续脚本需执行权限(x)才能运行）。
sudo chmod +x krakentools/bin/combine_kreports.py
# Kraken 报告转 MPA 格式：多个文件的Kraken reports转化为mpa文件
# 将 Bracken 校正后的 Kraken 报告（sample.S.bracken.kreport）转换为 MPA（Metaphlan）格式，便于后续多样性分析。
python krakentools/bin/kreport2mpa.py --display-header -r sample.S.bracken.kreport -o sample.S.bracken.MPA.TXT
# 合并所有样品的物种注释信息
# 将多个样品的 MPA 格式文件（如不同跳虫个体的肠道微生物数据）合并为一个矩阵，便于比较样品间的群落差异。
# 每个样品的 MPA 文件包含该样品中各物种的丰度，combine_mpa.py通过匹配物种名称（如k__Bacteria|s__XXX），将所有样品的丰度整合到一个矩阵中（行 = 物种，列 = 样品，单元格 = 丰度），为后续多样性分析提供输入。
python krakentools/bin/combine_mpa.py -i *.S.bracken.MPA.TXT -o combined.S.MPA.txt
# 获取alpha多样性指数
# 计算单个样品的 Alpha 多样性指数（此处为 Shannon 指数）。
# -f：输入 Bracken 校正后的丰度文件
python krakentools/bin/alpha_diversity.py -f sample.S.bracken -a Sh
# 获取beta多样性指数（Bray距离）
# 计算多个样品间的 Beta 多样性（此处为 Bray-Curtis 距离），反映样品间群落组成的差异。
# Beta 多样性用于衡量不同样品间的群落差异，Bray-Curtis 距离是常用指标，计算公式为：BC = 1 - [2×Σmin(pi, qi)] / (Σpi + Σqi)，其中pi和qi分别是两个样品中某物种的丰度。距离值范围为 0（完全相同）到 1（完全不同），值越大说明样品间群落差异越大。
python krakentools/bin/beta_diversity.py -i *.S.bracken --type bracken




# 使用ARGs_OAP进行ARGs比对（环境样本）
# 笔记：这里的ARGs是直接根据clean reads进行的比对。对环境样本进行了优化。灵敏度较高。支持read级别的分析。
# stage_one
args_oap stage_one -i input_reads_folder -o output_stage_one -f fastq -t 52
# stage_two
args_oap stage_two -i output_stage_one -o output_stage_two -t 52


# 使用rgi进行ARGs比对（肠道微生物）
# 需要先跑组装得到contigs，之后才能够跑rgi。因此后面在跑。








# HUMANN3功能比对（双端合成一个文件跑）
# 如果是双端测序，需要把双端文件CAT到一起再运行
cat sample_R1.fastq sample_R2.fastq > combined_sample.fastq
# 开始运行
humann -i combined_sample.fastq -o humann_output --threads 64
# 归一化丰度
humann_renorm_table -i humann_output/genefamilies.tsv --output humann_output/genefamilies_relab.tsv --units relab
humann_renorm_table -i humann_output/genefamilies.tsv --output humann_output/genefamilies_cpm.tsv --units cpm
# 合并多样品的表格到一起
humann_join_tables --input humann_output/genefamilies/normal --output humann_output/humann_2_genefamilies.tsv --file_name genefamilies
humann_join_tables --input humann_output/genefamilies/normal --output humann_output/humann_4_pathabundance.tsv --file_name pathabundance


# # HUMANN3功能比对（双端分开跑）
# mamba安装conda环境(conda渠道有问题)
conda activate base
mamba create -n humann3 -c conda-forge -c bioconda humann=3.6 python=3.10 -y
# 下载数据库
humann_databases
# 会显示最新的可以下载的地址。然后cd到目录，用wget下载。或者浏览器下载后ftp传过去。
# HUMANN3功能比对
# 将 KneadData 处理后的 clean reads 作为 HUMAnN3 的输入。对于双端测序数据，通常的做法是将 KneadData 输出的 *_paired_1.fastq.gz 和 *_paired_2.fastq.gz 文件合并成一个文件作为输入。虽然有些工作流程会自动处理双端，但合并输入是一个兼容性较好的方法。保留的 unpaired reads 也可以单独作为一个文件输入。
# 假设你的双端 clean reads 为 sample.R1.qualified.fastq.gz 和 sample.R2.qualified.fastq.gz
humann \
  --input sample.R1.qualified.fastq.gz,sample.R2.qualified.fastq.gz \
  --output humann_output \
  --threads 64

# 归一化：相对丰度（relative abundance, relab）
humann_renorm_table \
  --input humann_output/sample_genefamilies.tsv \
  --output humann_output/sample_genefamilies_relab.tsv \
  --units relab

# 归一化：每百万计数（CPM）
humann_renorm_table \
  --input humann_output/sample_genefamilies.tsv \
  --output humann_output/sample_genefamilies_cpm.tsv \
  --units cpm

# 合并多样品的表格到一起
humann_join_tables \
  --input humann_output/genefamilies/normal \
  --output humann_output/humann_2_genefamilies.tsv \
  --file_name genefamilies

humann_join_tables \
  --input humann_output/genefamilies/normal \
  --output humann_output/humann_4_pathabundance.tsv \
  --file_name pathabundance








######组装

# 组装（得到contigs）
# 合并不同平行的reads，用于混合组装
# cat sample1_R1.fastq sample2_R1.fastq sample3_R1.fastq > all_samples_R1.fastq
# cat sample1_R2.fastq sample2_R2.fastq sample3_R2.fastq > all_samples_R2.fastq

# Spades组装（准确度更高，适合肠道微生物）
# spades.py --meta -t 52 -m 500 --pe1-1 all_samples_R1.fastq --pe1-2 all_samples_R2.fastq -o spades_output
# 在宏基因组研究中（如你的肠道微生物项目），应始终使用 spades.py --meta，而不是默认的 SPAdes 模式。

# megahit组装
# megahit -t 32 -1 all_samples_R1.fastq -2 all_samples_R2.fastq -o megahit_output

# 用 CoverM 对每个样本定量 contig 覆盖度
# coverm contig -1 sample1_R1.fq.gz -2 sample1_R2.fq.gz -r coassembly/contigs.fasta --min-covered-fraction 0.75 -o abundances.tsv


######组装(全流程)######

# 第1步：合并 clean reads 用于 co-assembly
# 01_merge_reads_for_coassembly.sh

# 合并所有样本的 clean reads，仅用于后续 co-assembly
mkdir -p results/merged

cat clean/*_clean_R1.fastq.gz > results/merged/all.R1.fq.gz
cat clean/*_clean_R2.fastq.gz > results/merged/all.R2.fq.gz

echo "✅ Reads merged for co-assembly"

# # 第2步：metaSPAdes 全局组装
# # 02_spades_coassembly.sh
# # 需要大量内存mem，实测组内的服务器1T跑不了600G的合并reads。如果这个不行就换megahit。

# # 激活 SPAdes 环境
# source activate env_spades

# # 创建输出目录
# mkdir -p results/spades

# # 执行 metaSPAdes 组装
# spades.py \
#   --meta \
#   -1 results/merged/all.R1.fq.gz \
#   -2 results/merged/all.R2.fq.gz \
#   -t 48 \
#   -m 600 \
#   -o results/spades

# # 检查组装是否成功
# if [ ! -f results/spades/contigs.fasta ]; then
#   echo "❌ SPAdes failed: contigs.fasta not found"
#   exit 1
# fi

# echo "✅ Co-assembly completed"


# 第2步：megahit 全局组装
# 03_megahit_assembly.sh

# 创建输出目录
mkdir -p ./6_assembly/2_megahit

# 运行 MEGAHIT（专为大内存限制优化）
megahit \
  -1 ./6_assembly/1_mergereads/all.R1.fastq \
  -2 ./6_assembly/1_mergereads/all.R2.fastq \
  --presets meta-sensitive \
  --mem-flag 1 \
  -t 72 \
  --out-dir ./6_assembly/2_megahit \
  --out-prefix mg_6 \
  --keep-tmp-files \
  --tmp-dir ./6_assembly/2_megahit/tmp

# 验证结果
if [ ! -f ./6_assembly/2_megahit/fc_gut.contigs.fa ]; then
  echo "❌ MEGAHIT 组装失败！"
  exit 1
fi

echo "✅ MEGAHIT 完成！结果: ./6_assembly/2_megahit/fc_gut.contigs.fa"


# 第3步：过滤 ≥1 kb contigs（满足 OPERA-MS 要求）
# 03_filter_contigs_1k.sh

# 激活 seqkit 环境
source activate env_seqkit

# 创建输出目录
mkdir -p results/filtered

# 保留长度 ≥1000 bp 的 contigs
seqkit seq \
  -m 1000 \
  results/spades/contigs.fasta \
  > results/filtered/contigs_1k.fasta

echo "✅ Contigs ≥1 kb filtered"

# 第4步：估计文库插入片段大小
# 04_estimate_insert_size.sh
# bbmap不是conda程序,而是一堆sh脚本.放在了/database/work/zryan/software/bbTools/bbmap路径下.

# 创建临时目录
mkdir -p results/insert_size

# 提取前 10 万 reads 子集
zcat results/merged/all.R1.fq.gz | head -400000 | gzip > results/insert_size/subset.R1.fq.gz
zcat results/merged/all.R2.fq.gz | head -400000 | gzip > results/insert_size/subset.R2.fq.gz

# 将子集比对到 contigs
bbmap.sh \
  in=results/insert_size/subset.R1.fq.gz \
  in2=results/insert_size/subset.R2.fq.gz \
  ref=results/filtered/contigs_1k.fasta \
  out=results/insert_size/aligned.sam \
  t=12 \
  nodisk

# 计算 insert size 均值与标准差
python3 -c "
import pysam, numpy as np
sizes = []
with pysam.AlignmentFile('results/insert_size/aligned.sam', 'r') as sam:
    for read in sam:
        if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and read.reference_name == read.next_reference_name:
            isize = abs(read.template_length)
            if 100 < isize < 1000:
                sizes.append(isize)
mean = int(np.mean(sizes)) if sizes else 400
std = max(50, int(np.std(sizes))) if sizes else 100
with open('results/insert_size/insert_size.txt', 'w') as f:
    f.write(f'MEAN_INSERT={mean}\nSTD_INSERT={std}\n')
"

# 清理临时文件
rm results/insert_size/subset.R1.fq.gz results/insert_size/subset.R2.fq.gz results/insert_size/aligned.sam

echo "✅ Insert size estimated"

# 第5步：OPERA-MS scaffolding
# 05_opera_ms_scaffolding.sh

# 激活 OPERA-MS 环境
source activate env_opera_ms

# 加载 insert size 参数
source results/insert_size/insert_size.txt

# 创建输出目录
mkdir -p results/opera_ms

# 执行 scaffolding
opera_ms \
  --contig-file results/filtered/contigs_1k.fasta \
  --short-read1 results/merged/all.R1.fq.gz \
  --short-read2 results/merged/all.R2.fq.gz \
  --insert-size-mean $MEAN_INSERT \
  --insert-size-stdev $STD_INSERT \
  --num-processes 24 \
  --output-dir results/opera_ms

# 检查输出
if [ ! -f results/opera_ms/scaffolds.fasta ]; then
  echo "❌ OPERA-MS failed: scaffolds.fasta not found"
  exit 1
fi

echo "✅ Scaffolding completed"

# 第6步：QUAST 组装质量评估
# 06_quast_evaluation.sh

# 激活 QUAST 环境
source activate env_quast

# 创建输出目录
mkdir -p results/quast

# 比较 contigs 与 scaffolds 质量
quast.py \
  results/filtered/contigs_1k.fasta \
  results/opera_ms/scaffolds.fasta \
  -o results/quast \
  --threads 24 \
  --min-contig 1000 \
  --fast \
  --labels "Contigs_1k,Scaffolds"

echo "✅ QUAST report generated"

# 第7步：ORF 预测（在 scaffolds 上）
# 07_orf_prediction.sh

# 激活 Prodigal 环境
source activate env_prodigal

# 创建输出目录
mkdir -p results/orf

# 预测开放阅读框（meta 模式）
pprodigal \
  -i results/opera_ms/scaffolds.fasta \
  -o results/orf/genes.gff \
  -d results/orf/genes.fasta \
  -a results/orf/proteins.faa \
  -p meta \
  -f gff

echo "✅ ORF prediction completed on scaffolds"

# 第8步：构建非冗余基因集
# 08_build_nonredundant_gene_catalog.sh

# 激活 CD-HIT 环境
source activate env_cdhit

# 创建输出目录
mkdir -p results/gene_catalog

# 注意：当前流程仅有一个 scaffolds.fasta，故 genes.fasta 已是单一集合
# 此步骤主要用于兼容多组装合并场景，或未来扩展
cp results/orf/genes.fasta results/gene_catalog/combined_genes.fasta

# 使用 CD-HIT-EST 去冗余（90% identity）
cd-hit-est \
  -i results/gene_catalog/combined_genes.fasta \
  -o results/gene_catalog/nonredundant_genes.fasta \
  -M 0 \
  -c 0.9 \
  -n 8 \
  -T 160

echo "✅ Non-redundant gene catalog built"

# 第9步：KEGG 功能注释
# 09_kegg_annotation.sh

# 激活 KofamScan 环境
source activate env_kofamscan

# 创建输出目录
mkdir -p results/kegg ko_tmp

# 对蛋白质序列进行 KEGG Orthology 注释( -p -k 请替换为实际路径)
exec_annotation \
  results/orf/proteins.faa \
  -f mapper \
  -p /path/to/kofam_profiles \
  -k /path/to/kofam_ko_list \
  --cpu 60 \
  --tmp-dir ko_tmp \
  -o results/kegg/ko_output.txt

# 清理临时目录
rm -rf ko_tmp

echo "✅ KEGG annotation completed"

# 第10步：CoverM 样本级丰度定量
# 10_coverm_per_sample_quantification.sh

# 激活 CoverM 环境
source activate env_coverm

# 创建输出目录
mkdir -p results/coverm

# 获取样本列表（从 clean/ 目录）
samples=$(ls clean/*_clean_R1.fastq.gz | xargs -n1 basename | sed 's/_clean_R1.fastq.gz//')

# 对每个样本独立定量 scaffolds 覆盖度
for sample in $samples; do
  echo "Processing $sample..."
  coverm contig \
    --method metabat \
    --reference results/opera_ms/scaffolds.fasta \
    --forward clean/${sample}_clean_R1.fastq.gz \
    --reverse clean/${sample}_clean_R2.fastq.gz \
    --threads 8 \
    --min-covered-fraction 0.75 \
    --min-read-percent-identity 0.95 \
    --output-format sparse \
    > results/coverm/${sample}.tsv
done

# 合并为丰度矩阵
header="#scaffold"
for sample in $samples; do
  header="$header\t$sample"
done
echo -e "$header" > results/coverm/abundance_matrix.tsv

paste \
  <(cut -f1 results/coverm/${samples[0]}.tsv) \
  <(for s in $samples; do cut -f2 results/coverm/$s.tsv; done | paste -s -d$'\t') \
  >> results/coverm/abundance_matrix.tsv

echo "✅ Per-sample scaffold abundance matrix generated"

#####end#####




####原脚本补充内容
######预测ORF
# prodigal基因预测（预测出来ORF）
# pprodigal -i megahit_output/final.contigs.fa -o predicted_genes.gff -d predicted_genes.fasta -a predicted_proteins.faa -p meta -f gff


######KEGG功能基因比对
# 使用kofamscan进行比对
# exec_annotation predicted_proteins.faa -f mapper -p kofamscan_profiles -k kofamscan_ko_list --cpu 60 --tmp-dir ko_tmp -o ko_output.txt
# 这个predicted_proteins.faa文件是prodigal预测出来的蛋白质序列文件
# 根据KEGG官网上整理的KEGG不同level以及不同通路map进行vlookup函数比对

# ##### 非冗余基因集构建
# # 非冗余基因集构建
# # 把所有要分析的prodigal后contig的基因序列文件合并在一起
# cat predicted_genes.fasta other_sample_genes.fasta > combined_genes.fasta
# # CD-HIT去冗余(mmseqs2)
# cd-hit-est -i combined_genes.fasta -o nonredundant_genes.fasta -M 0 -c 0.9 -n 8 -T 160


##### 基因丰度比对(基于reads) #####
# 注意:非冗余基因集需要contigs跑到cd-hit之后才能够出来.因此虽然基于reads,但是需要在contig之后才能进行.
# 两个比对流程二选一,都是基于reads水平进行的丰度比对.

# 不基于比对的流程：使用salmon定量
# 构建索引
salmon index -t nonredundant_genes.fasta -i salmon_index
# 基因定量
salmon quant -i salmon_index -l A -p 60 --meta -1 sample_R1.fastq -2 sample_R2.fastq -o salmon_quant
# 生成基因count和TPM文件
salmon quantmerge --quants salmon_quant --column numreads -o gene_counts.txt
salmon quantmerge --quants salmon_quant --column tpm -o gene_tpm.txt
# 使用diamond获取基因注释表
# 构建数据库
diamond makedb --in nonredundant_genes.fasta -d nonredundant_genes_db
# 序列比对
diamond blastp -d nonredundant_genes_db -q query_proteins.fa -o blastp_output.txt --outfmt 6 --evalue 1e-5 --max-target-seqs 1 --threads 12

# 基于比对的流程：使用bwa定量
# 构建索引
bwa index -p bwa_index nonredundant_genes.fasta
# 比对
bwa mem -t 8 bwa_index sample_R1.fastq sample_R2.fastq -o bwa_output.sam
# 转换sam文件为bam文件并排序
samtools view -@ 20 -S -b -o bwa_output.bam bwa_output.sam
samtools sort -O bam -@ 40 -o bwa_output_sorted.bam bwa_output.bam
# 使用Bedtools进行丰度统计
bedtools bamtobed -i bwa_output_sorted.bam > bwa_output_sorted.bed
# 使用diamond获取基因注释表
# 构建数据库
diamond makedb --in nonredundant_genes.fasta -d nonredundant_genes_db
# 序列比对
diamond blastp -d nonredundant_genes_db -q query_proteins.fa -o blastp_output.txt --outfmt 6 --evalue 1e-5 --max-target-seqs 1 --threads 12

# 根据获取到的基因注释表，自行在excel中vlookup得到对应基因的丰度




# 定量结合salmon进行丰度计算

# diamond+card/vfdb/mobileOG-db
# diamond参数：
# identity
# coverage
# evalue

######分箱

# 使用MetaBAT2分箱
# 比对与排序
bowtie2-build --threads 72 assembled_contigs.fa assembled_contigs
bowtie2 -x assembled_contigs -1 sample_R1.fastq -2 sample_R2.fastq -p 20 -S alignment.sam
samtools view -@ 40 -b -S alignment.sam -o alignment.bam
samtools sort -@ 40 -l 9 -O BAM alignment.bam -o sorted_alignment.bam
jgi_summarize_bam_contig_depths --outputDepth contig_depth.txt sorted_alignment.bam
metabat2 -t 60 -i assembled_contigs.fa -a contig_depth.txt -o bins_dir/bin


######CoverM计算Contig，MAG等丰度
# 计算MAG丰度
coverm genome -d mag_directory -x fa -t 5 -c clean_reads_directory/*.fastq > mag_abundance_output --methods tpm

######GTDB-Tk对MAG进行物种比对及建树
# 物种注释
gtdbtk classify_wf --genome_dir mag_directory --out_dir classify_output --extension fa --mash_db mashdb --prefix mag --cpus 20
# 建树
gtdbtk de_novo_wf --genome_dir mag_directory --bacteria --outgroup_taxon p__Patescibacteria --out_dir tree_output --cpus 20

#######MAG携带的MGE比对

# 使用mobileOG-db比对
cd mag_directory
./mobileOGs-pl-kyanite.sh -i mag.fasta -d mobileOG-db/mobileOG-db_beatrix-1.6.All.dmnd -m mobileOG-db/mobileOG-db-beatrix-1.X.All.csv -k 1 -e 1e-20 -p 90 -q 90

######MAG的ARGs比对

# 使用diamond比对
diamond blastp \
  --threads 120 \
  -d sarg_database.dmnd \
  -q mag_proteins.fasta \
  -e 1e-10 \
  --id 90 \
  --query-cover 90 \
  --max-target-seqs 1 \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  -o arg_output.txt