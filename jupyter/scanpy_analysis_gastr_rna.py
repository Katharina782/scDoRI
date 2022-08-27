#!/usr/bin/env python
# coding: utf-8

# In[4]:


import scanpy as sc


# In[2]:


import numpy as np
import pandas as pd


# Note that rows correspond to cells and columns to genes. Eeach call to the `AnnData` object adds annotation to the expression matrix. `adata.obs` stores the annotation of observations and `adata.var` as `pd.DataFrame`. Unstructured annotations are saved as dictionaries in `adata.uns`. 

# In[3]:


adata = sc.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/anndata_rna.h5ad")


# In[4]:


results_file = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/scanpy_write_anndata_rna.h5ad"


# In[1]:


adata.X


# In[6]:


adata.obs["celltype.mapped"]


# In[7]:


adata.obs


# In[8]:


adata.uns["celltype.mapped_colors"]


# In[9]:


adata.var


# In[10]:


adata.X


# Normalize the data to 10, 000 reads per cell. This way counts become comparable between cells. Keepin mind that there might be dropout (zero counts due to sampling)

# ### Quality Control & Normalization
# 
# Quality control has already been performed for us, so we can just filter the object based on the "pass_rnaQC" column. Then we will normalize and scale the counts. 

# In[11]:


adata.obs["pass_rnaQC"] 


# It turns out that only cells which pass QC are included in the dataset, so we do not have to subset anymore.

# In[12]:


adata.obs[adata.obs["pass_rnaQC"] == True]


# In[13]:


# after normalization we want every cell to hae the same total counts
# therefore, we normalize by total counts over all genes
sc.pp.normalize_total(adata, target_sum=1e4)


# Logarithmize the data, because we are interested in relative changes.

# In[14]:


sc.pp.log1p(adata)


# In[15]:


adata


# Identify highly variable genes

# In[16]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[17]:


sc.pl.highly_variable_genes(adata)


# We set the raw attribute to the normalized and logarithmized data. This way we can go back later on. 

# In[18]:


adata.raw = adata


# In[ ]:





# In[19]:


# The highly variable genes are saved as follows:
adata.var.highly_variable


# Q: Should I use regress_out?

# In[20]:


# regress out effects of total counts per cell and mitochondrial gens
#sc.pp.regress_out(adata, ['nCount_RNA', 'mitochondrial_percent_RNA'])

# scale to unit variance
sc.pp.scale(adata, max_value=10)


# In[21]:


sc.tl.pca(adata, svd_solver='arpack')
    


# In[22]:


sc.pl.pca(adata, color = "Rp1")


# What is the contribution of single PCs to the total variance in the data? How many PCs should we consider for neighborhood graphs and clustering?

# In[23]:


sc.pl.pca_variance_ratio(adata, log=True)


# #### Compute the neighbourhood graph
# Which values should be used here? Try to play around with this

# Tried using 20 PCs for neighbor search.

# In[24]:


help(sc.pp.neighbors)


# In[25]:


sc.pp.neighbors(adata, n_pcs=20, knn=True, n_neighbors=15, method="umap")#, n_neighbors=10, n_pcs=40)


# In[26]:


adata.write(results_file)


# ### Clustering the neighborhood graph

# In[27]:


help(sc.tl.leiden)


# In[28]:


sc.tl.leiden(adata, resolution=.9)


# In[30]:


sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')


# In[32]:


sc.set_figure_params(figsize=(10,10))


# In[33]:


adata.obs["celltype.mapped"].unique()


# In[34]:


colPalette_celltypes = ['#532C8A',
 '#c19f70',
 '#f9decf',
 '#c9a997',
 '#B51D8D',
 '#3F84AA',
 '#9e6762',
 '#354E23',
 '#F397C0',
 '#ff891c',
 '#635547',
 '#C72228',
 '#f79083',
 '#EF4E22',
 '#989898',
 '#7F6874',
 '#8870ad',
 '#647a4f',
 '#EF5A9D',
 '#FBBE92',
 '#139992',
 '#cc7818',
 '#DFCDE4',
 '#8EC792',
 '#C594BF',
 '#C3C388',
 '#0F4A9C',
 '#FACB12',
 '#8DB5CE',
 '#1A1A1A',
 '#C9EBFB',
 '#DABE99',
 '#65A83E',
 '#005579',
 '#CDE088',
 '#f7f79e',
 '#F6BFCB']


# In[35]:


sc.pl.umap(adata, color=["leiden", "celltype.mapped", "nCount_RNA", "mitochondrial_percent_RNA", "sample"],
           title=["Leiden Clustering", "Mapped Celltypes", "nCount_RNA", "Percentage of mitochondrial genes", "Sample"], 
           ncols = 2, size=.9)


# In[36]:


#sc.tl.rank_genes_groups(adata, "leiden", method="logreg")
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[37]:


adata.write(results_file)


# In[38]:


adata.write_csvs(dirname="/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/anndata_as_csvs/", skip_data=False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#import anndata2ri as ar


# In[ ]:


from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter

