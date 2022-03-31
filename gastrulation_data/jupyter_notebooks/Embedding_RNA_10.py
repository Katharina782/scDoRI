#!/usr/bin/env python
# coding: utf-8

# ## Model with 10 dimensions in latent space

# In[1]:


import scanpy as sc
import scvi
import numpy as np


# In[2]:


adata = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/anndata_rna.h5ad")


# In[3]:


sc.pp.filter_genes(adata, min_counts=3)


# In[4]:


adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`
adata.raw


# In[5]:


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="sample"
)


# In[6]:


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["sample"],
    continuous_covariate_keys=["mitochondrial_percent_RNA", "ribosomal_percent_RNA"]
)


# In[7]:


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


# In[8]:


#model.train()


# In[ ]:





# In[9]:


# load anndata with latent space
adata = scvi.data.read_h5ad("gpu_trained_10_dim/anndata_object")


# In[10]:


adata


# In[12]:


# load trained model
model = scvi.model.SCVI.load("gpu_trained_10_dim/", adata=adata, use_gpu=False)


# In[ ]:


model


# In[13]:


latent = model.get_latent_representation()


# In[22]:


latent.shape


# In[ ]:


adata


# ## Visualization without batch correction

# In[14]:


sc.set_figure_params(figsize=(8,8))


# In[15]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.1)


# In[16]:


sc.pl.umap(
    adata,
    color=["celltype"],
    frameon=False,
    size=15
)
sc.pl.umap(
    adata,
    color=["sample"],
    ncols=2,
    frameon=False,
    size=15
)


# ## Visualization with batch correction 

# In[17]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)


# In[18]:


sc.pl.umap(
    adata,
    color=["celltype"],
    frameon=False,
    size=15
)
sc.pl.umap(
    adata,
    color=["sample"],
    ncols=2,
    frameon=False,
    size=15
)


# ## Have a look at the nearest neighbors embedding

# In[ ]:


adata.uns["neighbors"]


# In[ ]:


connections = adata.obsp["connectivities"]


# In[ ]:


connections = connections[0:500, 1:500].todense()

