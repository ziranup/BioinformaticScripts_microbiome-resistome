    # 运行 Prodigal
    # -i: 输入文件
    # -o: GFF 输出 (坐标信息)
    # -d: 核苷酸序列输出 (CDS)
    # -a: 氨基酸序列输出 (蛋白)
    # -p meta: 宏基因组模式 (适用于环境样本或单细胞)
    # -f gff: 输出格式
    # -q: 静默模式

# 假设 SARG.fasta 是多条完整 CDS（每条一个 ARG）
prodigal \
    -i /database/work/zryan/biodb/SARG/SARG_3.2.4/SARG.fasta \
    -o ./SARG_genes.gff \
    -d ./SARG_genes.fna \
    -a ./SARG_proteins.faa \
    -p meta \
    -q