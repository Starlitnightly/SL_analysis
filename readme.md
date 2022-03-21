# A pan-cancer atlas of tumour development of synthetic lethal gene developmental trajectory changes in tumour-infiltrating bone marrow cells 

This repository contains the scripts and function to reproduce the results of SL analysis

## Content

- `SL_Function.py`: the analysis function of SL
- `data`: the SL pair
- `script`: the pipeline of SL analysis

## Tutorial

A simple tutorial on SL analysis will be given here

```python
from SL_Function import SL_Analysis
import anndata
#import data
adata=anndata.read_h5ad("data/lym_dpt.h5ad")
SL_pd=pd.read_csv('data/sl-for-scrnaseq-full.csv')
#Analysis Cancer 
adata_c=adata[adata.obs['tissue']=='T']
SL_c=SL_Analysis(adata_c,SL_pd)
sl_c.Lazy_analysis(['SL_count.csv'.format(cancer,trajectory),'SL_slope.csv'.format(cancer,trajectory),'SL_pair.csv'.format(cancer,trajectory)])
#Analysis Paracancer
adata_p=adata[adata.obs['tissue']=='N']
SL_p=SL_Analysis(adata_p,SL_pd)
sl_p.Lazy_analysis(['SL_count.csv'.format(cancer,trajectory),'SL_slope.csv'.format(cancer,trajectory),'SL_pair.csv'.format(cancer,trajectory)])
```

## Data

The parsed data can be downloaded in each folder.

## Contact

- Zehua Zeng ([starlitnightly@163.com](mailto:starlitnightly@163.com))