import scanpy as sc
import numpy as np

from .. import logging as logg
from .. import settings

import shutil

try:
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    import rpy2.robjects as robjects

    pandas2ri.activate()
except Exception as e:
    raise Exception(
        "rpy2 installation is necessary for "
        + task
        + '. \
        \nPlease use "pip3 install rpy2" to install rpy2'
    )

if not shutil.which("R"):
    raise Exception(
        "R installation is necessary."
    )

try:
    rstats = importr("stats")
except Exception as e:
    raise Exception(
        "R installation is necessary."
    )

try:
    SE = importr("SummarizedExperiment")
except Exception as e:
    raise Exception(
        'R package "SummarizedExperiment" is necessary.\n'
        'Please install it from https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html and try again'
    )    
    
    
try:
    deseq = importr("DESeq2")
except Exception as e:
    raise Exception(
        'R package "DESeq2" is necessary.\n'
        'Please install it from https://bioconductor.org/packages/release/bioc/html/DESeq2.html and try again'
    )
    
def run(adata,formula,**kwargs):
    logg.info("Running DESeq2", reset=True, end="\n")
    
    adata.uns["Formula"] = formula 
    if isinstance(adata.X, np.ndarray) == False:
        logg.warn("Sparse matrix detected, densifying...")
        adata.X = adata.X.A
    dds = deseq.DESeqDataSetFromMatrix(countData=adata.X.T, 
                                        colData=adata.obs,
                                        design=Formula(formula))
    dds = deseq.DESeq(dds, **kwargs)
    adata.layers["normalized"] = deseq.counts_DESeqDataSet(dds, normalized=True).T 
    
    adata.uns["dds"] = dds
    
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .layers['normalized'] normalized count matrix.\n"
        "    .uns['dds'] DESeq2 R object.\n"
        "    .uns['Formula'] formula used for design parameter."
    )
    
    
def vst(adata):
    logg.info("Obtaining vsd", end="\n")
    vsd=deseq.vst(adata.uns["dds"], blind=False)
    adata.layers["vsd"] = SE.assay(vsd).T
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .layers['vsd'] variance stabilized count matrix."
    )
    
def pca(adata,top_genes=500):
    logg.info("Obtaining PCA", reset=True, end="\n")
    raw = adata.X.copy() 
    mat = "vsd" if "vsd" in adata.layers else "normalized"
    adata.X = adata.layers[mat]
    adata.var[mat+"_std"]=adata.X.std(axis=0)
    adata.var["highly_variable"]=False
    adata.var.loc[adata.var[mat+"_std"].sort_values(ascending=False).index[:top_genes],
                  "highly_variable"]=True
    logg.info("    on highly variable genes using "+mat+ " matrix", end="\n")
    sc.pp.pca(adata)
    adata.X = raw
    
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .var['"+mat+"_std'] variance of genes calculated from vsd matrix.\n"
        "    .var['highly_variable'] genes considered as highly variable.\n"
        "    .obsm['X_pca'] PCA results.\n"
        "    .uns['pca'] PCA additional results."
    )

def show_results(adata):
    return list(deseq.resultsNames(adata.uns["dds"]))

def result(adata,name,lfc_shrink=False,**kwargs):
    
    if ~np.isin(name,show_results(adata)):
        raise ValueError('result not present in DESeq2 object, available results are:\n    '+
                         '\n    '.join(show_results(adata)))
    logg.info("Generating DE results", reset=True, end="\n")
    to_dataframe = robjects.r('function(x) data.frame(x)')
    if lfc_shrink:
        logg.info("    running LFC shrinking", end="\n")
        df=to_dataframe(deseq.lfcShrink(adata.uns["dds"],coef=name))
        df=robjects.conversion.rpy2py(df)
        df.index=adata.var_names
        if name in adata.uns:
            adata.uns[name]["LFC_shrink"]=df
        else:
            adata.uns[name]={"LFC_shrink":df}
    else:
        df=to_dataframe(deseq.results(adata.uns["dds"],name=name))
        df=robjects.conversion.rpy2py(df)
        df.index=adata.var_names
        if name in adata.uns:
            adata.uns[name]["Results"]=df
        else:
            adata.uns[name]={"Results":df}
            
    res = "LFC_shrink" if lfc_shrink else "Results"
    logg.info(
        "    done",
        time=True,
        end=" " if settings.verbosity > 2 else "\n",
    )
    logg.hint(
        "added\n"
        "    .uns['"+name+"']['"+res+"'] table of differential expression results."
    )
            
    
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
        "    "+name+".h5ad: anndata file containing counts matrices, pca and DE results.\n"
        "    "+name+".rds: R object of the DESeqDataSet."
    )

    
def read(name):
    adata = sc.read(name+".h5ad")
    readRDS = robjects.r['readRDS']
    adata.uns["dds"] = readRDS(name+'.rds')
    return adata
