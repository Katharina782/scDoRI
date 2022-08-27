#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import scvi
import numpy as np
import anndata as ad
from anndata import AnnData
from scipy import sparse
import pandas as pd


# In[2]:


#multiome = scvi.data.read_h5ad("multiome_gastr.h5ad")


# In[3]:


rna_adata = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/anndata_rna.h5ad")
atac_adata = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/anndata_atac_peak_matrix.h5ad")


# In[4]:


#adata = scvi.data.read_10x_multiome("filtered_feature_bc_matrix")


# In[5]:


#adata.var_names_make_unique()
rna_adata.var_names_make_unique()
atac_adata.var_names_make_unique()


# In[6]:


old = scvi.data.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/old_data_three_timepoints/anndata_rna.h5ad")


# In[7]:


rna_adata = rna_adata[:,rna_adata.var.index.isin(old.var.index)]


# In[8]:


annotations = pd.merge(old.var, rna_adata.var, on="gene", how = "inner")
annotations.index = old.var.index


# In[9]:


annotations = annotations.drop(["gene", "Strand"], axis=1)


# In[10]:


annotations = annotations.rename(columns={"Accession":"ID", "End":"end", "Start":"start", "Chromosome":"chr"})


# In[11]:


rna_adata.var = annotations


# In[12]:


matrix = sparse.hstack((rna_adata.X, atac_adata.X))


# In[13]:


rna_df = rna_adata.var
#rna_df = rna_df.rename(columns = {"gene":"ID"})
rna_df["modality"] = ["Gene_expression" for i in range(len(rna_df))]
rna_df


# In[14]:


atac_df = atac_adata.var
atac_df = atac_df.rename(columns={"idx":"ID"})
atac_df["modality"]=["Peaks" for i in range(len(atac_df))]


# In[15]:


atac_df = atac_df.drop(["score"], axis = 1)


# In[16]:


atac_df


# In[17]:


df = rna_df.append(atac_df)


# In[18]:


multiome = AnnData(
X = matrix,
obs = atac_adata.obs,
var = df,
uns = rna_adata.uns)


# In[19]:


multiome.shape


# In[20]:


adata_mvi = scvi.data.organize_multiome_anndatas(multiome)


# In[21]:


adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()


# In[22]:


print(adata_mvi.shape)
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
print(adata_mvi.shape)


# In[23]:


scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key='modality')


# In[26]:


mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var['modality']=='Gene_expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
)
mvi.view_anndata_setup()


# In[ ]:


mvi.train()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#multiome_mvi = scvi.data.organize_multiome_anndatas(multiome)


# In[ ]:


#multiome_mvi.var


# In[ ]:


#adata_mvi = scvi.data.organize_multiome_anndatas(multiome)


# In[ ]:


#multiome.write("multiome_gastr.h5ad", compression="gzip")


# In[ ]:





# In[ ]:





# In[ ]:


#multiome_mvi = multiome_mvi[:, multiome_mvi.var["modality"].argsort()].copy()


# In[ ]:


#multiome_mvi.var


# Filter features which appear in fewer than 1% of cells.

# The main batch key should correspond to the modality of the cells. All other batch annotations are added via categorical_covariate_keys.

# In[ ]:


test = scvi.data.organize_multiome_anndatas(multiome)


# In[ ]:


scvi.model.MULTIVI.setup_anndata(adata)#,
                                # batch_key = "modality")
                                 #categorical_covariate_keys=["sample"])


# We need to specify how many features are genes and how many are genomic regions. 

# In[ ]:


mvi = scvi.model.MULTIVI(
    adata,
    n_genes=(adata.var['modality']=='Gene_expression').sum(),
    n_regions=(adata.var['modality']=='Peaks').sum()
    #dropout_rate=0.2,
    #n_layers_encoder=2,
    #n_layers_decoder=2,
    #n_latent=20
)


# arches_params = dict(
#     n_genes=(multiome.var["modality"]=="Gene_expression").sum(),
#     n_regions=(multiome.var["modality"]=="Peaks").sum(),
#     #dropout_rate=0.2,
#     #n_layers_encoder=2,
#     #n_layers_decoder=2,
#     #n_latent=20,
#     #latent_distribution='normal'
# )
# 
# mvi = scvi.model.MULTIVI(
#     multiome,
#     **arches_params
# )
# mvi.view_anndata_setup()

# In[ ]:


mvi.train()

