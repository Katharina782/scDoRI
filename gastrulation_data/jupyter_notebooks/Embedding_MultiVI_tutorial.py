#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import scvi
import numpy as np


# In[ ]:


import anndata as ad
from anndata import AnnData


# In[ ]:


from scipy import sparse


# In[ ]:


import pandas as pd


# In[50]:


#rna_adata = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/anndata_rna.h5ad")
#atac_adata = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/anndata_atac_peak_matrix.h5ad")


# In[ ]:


adata = scvi.data.read_10x_multiome("filtered_feature_bc_matrix")


# In[ ]:


adata


# In[ ]:





# In[ ]:


#multiome_mvi.var


# Filter features which appear in fewer than 1% of cells.

# The main batch key should correspond to the modality of the cells. All other batch annotations are added via categorical_covariate_keys.

# In[23]:


scvi.model.MULTIVI.setup_anndata(adata)#,
                                # batch_key = "modality")
                                 #categorical_covariate_keys=["sample"])


# We need to specify how many features are genes and how many are genomic regions. 

# In[29]:


mvi = scvi.model.MULTIVI(
    adata,
    n_genes=(adata.var['modality']=='Gene Expression').sum(),
    n_regions=(adata.var['modality']=='Peaks').sum(),
    n_layers_encoder=2,
    n_layers_decoder=2,
    n_latent=20,
)


# In[ ]:


mvi.train()

