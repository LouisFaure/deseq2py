import scanpy as sc
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from . import logging as logg
from . import settings

pandas2ri.activate()
deseq = importr('DESeq2')
SE = importr('SummarizedExperiment')


def run(adata,formula,top_genes=500,**kwargs):
    
    logg.info("Running DESeq2 pipeline", reset=True, end="\n")
    
    adata.uns["Formula"] = formula 
    dds = deseq.DESeqDataSetFromMatrix(countData=adata.X.T, 
                                        colData=adata.obs,
                                        design=Formula(formula))
    dds = deseq.DESeq(dds, **kwargs)
    adata.layers["normalized"] = deseq.counts_DESeqDataSet(dds, normalized=True).T
    
    logg.info("    obtaining vsd", end="\n")
    
    vsd=deseq.vst(dds, blind=False)
    raw = adata.X.copy()
    adata.X = SE.assay(vsd).T
    
    adata.var["vsd_std"]=adata.X.std(axis=0)
    adata.var["highly_variable"]=False
    adata.var.loc[adata.var["vsd_std"].sort_values(ascending=False).index[:top_genes],
                  "highly_variable"]=True
    
    logg.info("    PCA on highly variable genes", end="\n")
    
    sc.pp.pca(adata)
    
    adata.layers["vsd"] = adata.X
    adata.X = raw
    
    adata.uns["dds"] = dds
    
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .layers['normalized'] normalized count matrix.\n"
        "    .layers['vsd'] variance stabilized count matrix.\n"
        "    .obsm['X_pca'] PCA results.\n"
        "    .var['vsd_std'] variance of genes calculated from vsd matrix.\n"
        "    .var['highly_variable'] genes considered as highly variable.\n"
        "    .uns['dds'] DESeq2 R object.\n"
        "    .uns['pca'] PCA additional results.\n"
        "    .uns['Formula'] formula used for design parameter."
    )
    

def result(adata,name,lfc_shrink=False,**kwargs):
    logg.info("Generating DE results", reset=True, end="\n")
    to_dataframe = robjects.r('function(x) data.frame(x)')
    if lfc_shrink:
        logg.info("    running LFC shrinking", end="\n")
        df=to_dataframe(deseq.lfcShrink(adata.uns["dds"],coef=name))
        df.index=adata.var_names
        if "LFC_shrink" in adata.uns:
            adata.uns["LFC_shrink"][name]=df
        else:
            adata.uns["LFC_shrink"]={name:df}
    else:
        df=to_dataframe(deseq.results(adata.uns["dds"],name=name))
        df.index=adata.var_names
        if "Results" in adata.uns:
            adata.uns["Results"][name]=df
        else:
            adata.uns["Results"]={name:df}
            
    res = "LFC_shrink" if lfc_shrink else "Results"
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .uns['"+res+"']['"+name+"'] table of differential expression results."
    )
            
def plot(adata,mode,name):
    df=adata.uns[mode][name]
    fig,ax=plt.subplots(figsize=(4,3.5))
    ax.scatter(df.baseMean,df.log2FoldChange,s=2,c="lightgrey")
    ax.scatter(df.loc[df.padj<0.05].baseMean,df.loc[df.padj<0.05].log2FoldChange,s=2)
    ax.semilogx();
    ax.set_ylabel("log fold change");
    ax.set_xlabel("mean of normalized counts");
    ax.axhline(c="grey",linewidth=2)
    ax.grid(False)
    ax.set_title(name)
    
def save(adata,name):
    logg.info("Saving results", reset=True, end="\n")
    saveRDS = robjects.r['saveRDS']
    saveRDS(adata.uns["dds"], name+'.rds')
    del adata.uns["dds"]
    adata.write(name+'.h5ad')
    logg.info(
        "    done",
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "saved\n"
        "    "+name+".h5ad anndata file containing counts matrices, pca and DE results.\n"
        "    "+name+".rds R object of the DESeqDataSet."
    )