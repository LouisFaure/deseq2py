# DESeq2py
Python wrapper of DESeq2 R package, combined with anndata object for easy access and plotting of the data.

## Installation

The best way to get everything running smoothly is to use conda:

```bash
conda create -n deseq2py -c conda-forge -c bioconda python=3.8 bioconductor-deseq2 bioconductor-apeglm rpy2 -y
conda activate deseq2py
pip install git+https://github.com/LouisFaure/deseq2py.git
```
    
## Usage

```python
import scanpy as sc
import deseq2py as deseq
adata_bulk = sc.read("adata_bulk.h5ad")

deseq2.tl.run(adata_bulk,formula="~ condition")
deseq2.tl.vst(adata_bulk)
deseq2.tl.pca(adata_bulk)

deseq2.pl.pca(adata_bulk,color="condition",
              show_ellipses=True)

deseq2.tl.result(adata_bulk,lfc_shrink=True,name="condition_A_vs_B")
```
