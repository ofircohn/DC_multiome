# Figure 4 -scVelo
# load packages
import pandas as pd
from scipy.io import mmread
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import scanpy as sc
import anndata as ad
import pickle
import numpy as np
import scvelo as scv
import cellrank as cr
import multiprocessing
import matplotlib.colors
plt = get_matplotlib_pyplot(False, raise_if_not_available=True)
from matplotlib.colors import  LinearSegmentedColormap
coresNum = multiprocessing.cpu_count()


# Subset cells for velocity
adata = scv.read('/broad/sankaranlab/cohn/data/GRN/velocyto/scVelo_annData_DC.h5ad')
#subset = ['MEP / Ery / MAST_1', 'pDC_5']
subset = ['MEP / Ery / MAST_1', 'cDC1_1']
adata = adata[adata.obs['CellTypeAll'].isin(subset)]

# scVelo pipline
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs = 10) 
scv.tl.velocity(adata, mode = 'dynamical',n_jobs = 20)
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)
scv.tl.velocity_pseudotime(adata)

adata.obs.to_csv("cDC1_meta.csv")
