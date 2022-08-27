#!/usr/bin/env python
# coding: utf-8

# In[4]:


import scanpy as sc


# In[8]:


import numpy as np
import pandas as pd


# In[6]:


rna = sc.read_h5ad("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/anndata_rna.h5ad")


# In[15]:


matrix = rna.X.todense()


# In[16]:


np.int_(matrix)


# In[18]:


pd.DataFrame(np.int_(matrix)).to_csv("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/count_matrix_gastr.csv")


# In[21]:


import csv


# In[23]:


test = csv.reader("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/count_matrix_gastr.csv")


# In[ ]:


pd.read_csv("/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/count_matrix_gastr.csv")

