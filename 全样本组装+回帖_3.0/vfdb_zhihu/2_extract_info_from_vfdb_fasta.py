import re
import pandas as pd
# setA
with open('./data/SetA_anno.txt',"r") as f:
    data=f.readlines()
anno_list=[]
for line in data:
    info = re.findall('>(\S+)\(gb\|\S+\) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    if len(info)==0:
        info = re.findall('>(\S+) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    tmp = pd.DataFrame(info)
    tmp.columns=['VF_gene_id','VF_gene_symbol','Gene_description','VF_name','VF_id','VF_category_level1','VF_category_id','taxonomy']
    anno_list.append(tmp)
anno = pd.concat(anno_list)
f.close()
anno.to_csv('./SetA_info.tsv',sep='\t',index=0)    
#setB
import re
import pandas as pd
with open('./data/SetB_anno.txt',"r") as f:
    data=f.readlines()
anno_list=[]
for line in data:
    info = re.findall('>(\S+)\(gb\|\S+\) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    if len(info)==0:
        info = re.findall('>(\S+) \((.*?)\) (.*) \[(.*) \((.*)\) - (.*) \((.*)\)\] \[(.*)\]',line)
    tmp = pd.DataFrame(info)
    tmp.columns=['VF_gene_id','VF_gene_symbol','Gene_description','VF_name','VF_id','VF_category_level1','VF_category_id','taxonomy']
    anno_list.append(tmp)
anno = pd.concat(anno_list)
f.close()
anno.to_csv('./SetB_info.tsv',sep='\t',index=0) 