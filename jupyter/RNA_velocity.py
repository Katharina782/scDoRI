#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import scvelo as scv
#import cellrank as cr
import numpy as np
import pandas as pd
#import anndata as ad


# In[2]:


adata = adata = sc.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/anndata_rna.h5ad")


# In[3]:


adata


# In[4]:


help(scv.pp.moments)


# In[5]:


scv.pp.filter_and_normalize(adata)
# computes PCA and neighbors
scv.pp.moments(adata)


# Lets have a look at the proportions of spliced and unspliced RNA.

# In[6]:


scv.pl.proportions(adata)


# For each gene a steady state ratio of spliced and unspliced mRNA counts is fitted. Velocities are obtained as residuals from this ratio. The computed velocities are stored in `adata.layers`. Velocities are vectors in gene expression space.

# In[ ]:


scv.tl.velocity(adata, mode='stochastic')


# The combination of velocities across genes can be used to determine future states of an individual cell. To project the velocities into a lower dimensional space we compute cell-to-cell transition probabilities. For each velocity vectors the cell transtitions which fit to the direction are found. A transition probability is computed with cosine correlation between cell-to-cell transitions and velocity vectors. These cosine correlations are stored as a matrix = velocity graph $n_{obs} x n_{obs}$. The velocity graph summarizes the possible cell state changes and that are well explained through the velocity vectors. 

# In[ ]:


# project the velocities into lower dimensional space
scv.pl.velocity_embedding_stream(adata, basis='umap')

