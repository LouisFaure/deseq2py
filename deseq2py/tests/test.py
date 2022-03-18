import numpy as np
import pandas as pd
import scanpy as sc
import deseq2py as deseq2

cnts=pd.read_csv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152774/suppl/GSE152774_HTSeq_genes_count_df.csv.gz",index_col=0)
cnts.columns=cnts.columns.astype(str)

X=cnts.iloc[:,[0,2,9,1,4,6,3,5,7]].T.astype(np.float32)

adata=sc.AnnData(X)

adata.obs["conditions"]=["A","A","A","B","B","B","C","C","C"]
deseq2.tl.run(adata,formula="~ conditions")
deseq2.tl.vst(adata)
deseq2.tl.pca(adata)

deseq2.tl.show_results(adata)
deseq2.tl.result(adata,lfc_shrink=True,name="conditions_B_vs_A")

deseq2.pl.result(adata,mode="LFC_shrink",name="conditions_B_vs_A")
deseq2.pl.pca(adata,color="conditions", show_ellipses=True)

deseq2.tl.save(adata,"deseq2")
adata=deseq2.tl.read("deseq2")
