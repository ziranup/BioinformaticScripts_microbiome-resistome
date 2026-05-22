import pandas as pd 


VF_anno = pd.read_csv('./data/VFs.txt',encoding='unicode_escape',sep='\t')
VF_anno = VF_anno.loc[:,['VFID','VF_FullName','Structure','Function','Mechanism']]
VF_anno.columns=['VF_id','VF_FullName','Structure','Function','Mechanism']

# VFDB_Set_A_anno
alig = pd.read_csv('./diamond_SetA_out.tsv',sep='\t',header=None)
fasta_anno = pd.read_csv('./SetA_info.tsv',sep='\t')
alig = alig.iloc[:,[0,1]]
alig.columns = ['GENE','VF_gene_id']
alig['VF_gene_id'] = alig['VF_gene_id'].map(lambda x :x.split('(')[0])
VFDB_anno=pd.merge(alig,fasta_anno,on='VF_gene_id',how='inner')
VFDB_anno = pd.merge(VFDB_anno,VF_anno,on='VF_id',how='left')
VFDB_anno.to_csv('./VFDB_anno_SetA.tsv',sep='\t',index=0)

# VFDB_Set_B_anno
alig = pd.read_csv('./diamond_SetB_out.tsv',sep='\t',header=None)
fasta_anno = pd.read_csv('./SetB_info.tsv',sep='\t')
alig = alig.iloc[:,[0,1]]
alig.columns = ['GENE','VF_gene_id']
alig['VF_gene_id'] = alig['VF_gene_id'].map(lambda x :x.split('(')[0])
VFDB_anno=pd.merge(alig,fasta_anno,on='VF_gene_id',how='inner')
VFDB_anno = pd.merge(VFDB_anno,VF_anno,on='VF_id',how='left')
VFDB_anno.to_csv('./VFDB_anno_SetB.tsv',sep='\t',index=0)