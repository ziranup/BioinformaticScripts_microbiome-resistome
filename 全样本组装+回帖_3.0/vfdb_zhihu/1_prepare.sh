#!/bin/bash
set -e


# 下载
mkdir data
cd data
wget http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz

# 回到流程主目录
cd ..

# diamond建库与比对
# SetA
diamond makedb --in ./data/VFDB_setA_pro.fas.gz -d SetA 1>stdout_SetA 2>stderr_SetA
diamond blastp --db SetA -q /data/yanziran/mg_6_new/D_b2_CoORFs/proteins.faa -o diamond_SetA_out.tsv --max-target-seqs 1 --evalue 1e-7 --min-score 60 --block-size 40.0 --index-chunks 1 1>blastp_SetA.out 2>blastp_SetA.err --sensitive --threads 72
# SetB
diamond makedb --in ./data/VFDB_setB_pro.fas.gz -d SetB 1>stdout_SetB  2>stderr_SetB 
diamond blastp --db SetB  -q /data/yanziran/mg_6_new/D_b2_CoORFs/proteins.faa -o diamond_SetB_out.tsv --max-target-seqs 1 --evalue 1e-7 --min-score 60 --block-size 40.0 --index-chunks 1 1>blastp_SetB.out 2>blastp_SetB.err --sensitive --threads 72

# 从fasta中提取信息
zcat ./data/VFDB_setA_pro.fas.gz | grep '^>' > ./data/SetA_anno.txt
zcat ./data/VFDB_setB_pro.fas.gz | grep '^>' > ./data/SetB_anno.txt

# VFs.xls.gz转为文本
# 因为VFs.xls的格式不适合用python读入，会报错，因此我们将第一行删掉后，另存为tsv文件，然后读入，读入这里使用的是encoding='unicode_escape'，要不然也会读不进去，我也不知道为什么。这种编码的问题一直搞不太明白。