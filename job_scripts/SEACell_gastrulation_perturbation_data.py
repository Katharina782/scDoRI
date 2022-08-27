#!/usr/bin/env python
# coding: utf-8

# # Computing SEACells based on scATAC
# 
# Since in the notebooks and tutorials they recommend using ATAC data for computing metacells this is what I will do here. Instead of using SVD I will use the 10 dimensional latent embedding obtained with PeakVI.

# # Imports

# In[3]:


import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns




data_dir = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data_new/Kathi/"


# scRNA 
ad_rna =  sc.read_h5ad(data_dir + "rna_anndata_from_seurat")



# scATAC
ad_atac = sc.read_h5ad(data_dir + "perturbation_sea.h5ad")




ad_atac.layers["csr"] = ad_atac.X




#SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='csr')


atac_meta_ad, rna_meta_ad = SEACells.genescores.prepare_multiome_anndata(ad_atac, ad_rna, SEACell_label='SEACell')


sc.write(data_dir + "atac_pert_sea_agg.h5ad", adata = atac_meta_ad)

sc.write(data_dir + "rna_pert_sea_agg.h5ad", adata = rna_meta_ad)

