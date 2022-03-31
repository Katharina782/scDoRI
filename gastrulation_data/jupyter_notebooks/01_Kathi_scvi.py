#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc


# In[2]:


import scvi


# In[3]:


import numpy as np


# In[4]:


import pandas as pd


# # RNA

# In[21]:


adata = scvi.data.read_h5ad("anndata_rna.h5ad")


# Remove genes which are exprssed at very low levels or other outliers. Scvi tools use the raw count matrix. We can save the raw count matrix in an extra layer and then normalize the count matrix for other tasks. 

# In[6]:


sc.pp.filter_genes(adata, min_counts=3)


# In[7]:


adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`
adata.raw


# The counts matrix contains the raw counts:

# Reduce the number of features which will be used as input for the model. You can use 1,000 to 10,000 most highly variable genes. Lets try with 10,000 first. Here we can already add a batch_key which means that the highly variable genes will be selected for each "batch" separately. We will not select all highly variable genes from just one "batch". In our case the only batch effect we would expect is between cells of different replicates. 

# In[8]:


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="sample"
)


# In[9]:


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["sample"]
)


# In[10]:


arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    
    dropout_rate=0.2,
    n_layers=2,n_latent=20,gene_likelihood='nb'
)

model = scvi.model.SCVI(
    adata,
    **arches_params
)


# Create a model. 

# In[12]:


model


# Train the model.

# In[13]:


#model.train()


# Save the model

# In[14]:


#model.save("scvi_models/")


# In[15]:


model = scvi.model.SCVI.load("Manu_model/", adata=adata, use_gpu=False)


# ### scVI Output

# In[ ]:


latent = model.get_latent_representation()


# We can save the output of the model in the anndata object:

# In[ ]:


adata.obsm["X_scVI"] = latent


# # Visualization without batch correction

# In[ ]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.1)


# In[ ]:


sc.pl.umap(
    adata,
    color=["celltype"],
    frameon=False,
)
sc.pl.umap(
    adata,
    color=["sample"],
    ncols=2,
    frameon=False,
)


# # Visualization after Batch correction

# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.pl.umap(
    adata,
    color=["celltype"],
    frameon=False,
)
sc.pl.umap(
    adata,
    color=["sample"],
    ncols=2,
    frameon=False,
)

