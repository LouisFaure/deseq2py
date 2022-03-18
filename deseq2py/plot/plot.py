from typing import Optional, Union
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scanpy.plotting._utils import savefig_or_show
from .utils import confidence_ellipse

def pca(adata,
        color: Optional[str] = None,
        s: int = 200,
        show_ellipses: bool = False,
        ellipse_std: float = 3.0,
        ax: Optional = None,
        show: Optional[bool] = None,
        save: Union[str, bool, None] = None,
        **kwargs):
    
    
    if ax is None:
        ax = sc.pl.pca(adata,s=s,color=color,show=False,**kwargs)
    else:
        sc.pl.pca(adata,s=s,color=color,show=False,ax=ax,**kwargs)

    var=np.array([int(round(p)) for p in adata.uns["pca"]["variance_ratio"]*100])

    ax.set_ylabel("PC2 ({}% variance)".format(var[1]))
    ax.set_xlabel("PC1 ({}% variance)".format(var[0]))

    if show_ellipses:
        dct=dict(zip(adata.obs[color].cat.categories,adata.uns[color+"_colors"]))
        for c in adata.obs[color].cat.categories:
            pc=adata[adata.obs[color]==c].obsm["X_pca"]
            confidence_ellipse(pc[:,0],pc[:,1],ax,n_std=ellipse_std,
                               facecolor=dct[c],alpha=.1,
                               edgecolor=dct[c],zorder=-1)
            confidence_ellipse(pc[:,0],pc[:,1],ax,n_std=ellipse_std,
                               edgecolor=dct[c],zorder=-1)
            
    savefig_or_show("pca", show=show, save=save)

    if show == False:
        return ax

def result(adata,
           mode: str,
           name: str,
           figsize: tuple = (4,3.5),
           ax: Optional = None,
           show: Optional[bool] = None,
           save: Union[str, bool, None] = None):
    
    df=adata.uns[name][mode]
    if ax is None:
        fig,ax=plt.subplots(figsize=figsize,constrained_layout=True)
    ax.scatter(df.baseMean,df.log2FoldChange,
               s=2,c="lightgrey",rasterized=True)
    ax.scatter(df.loc[df.padj<0.05].baseMean,df.loc[df.padj<0.05].log2FoldChange,
               s=2,rasterized=True)
    ax.semilogx();
    ax.set_ylabel("log fold change");
    ax.set_xlabel("mean of normalized counts");
    ax.axhline(c="grey",linewidth=2)
    ax.grid(False)
    ax.set_title(name)
    
    savefig_or_show("result", show=show, save=save)

    if show == False:
        return ax