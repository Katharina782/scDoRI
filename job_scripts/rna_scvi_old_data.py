import scanpy as sc


# In[2]:


import scvi


# In[3]:


import numpy as np


# # RNA

# In[4]:

dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"
adata = scvi.data.read_h5ad(dir_data + "old_anndata_rna.h5ad")

adata = adata[adata.obs["celltype.mapped"].isin(["Forebrain_Midbrain_Hindbrain", "Neural_crest", "Rostral_neurectoderm", "Spinal_cord", "Caudal_neurectoderm", "Surface_ectoderm", "NMP", "Epiblast", "Caudal_epiblast"])] 


# Remove genes which are exprssed at very low levels or other outliers. Scvi tools use the raw count matrix. We can save the raw count matrix in an extra layer and then normalize the count matrix for other tasks. 

# In[5]:


sc.pp.filter_genes(adata, min_counts=3)


# In[6]:


adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`


# Reduce the number of features which will be used as input for the model. You can use 1,000 to 10,000 most highly variable genes. Lets try with 10,000 first. Here we can already add a batch_key which means that the highly variable genes will be selected for each "batch" separately. We will not select all highly variable genes from just one "batch". In our case the only batch effect we would expect is between cells of different replicates. 

# In[7]:


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="sample"
)


# In[12]:


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["sample"]
)


# In[13]:


arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    
    dropout_rate=0.2,
    n_layers=2,
    n_latent=10,
    gene_likelihood='nb'
)

model = scvi.model.SCVI(
    adata,
    **arches_params
)


# Create a model. 

# In[46]:


#model = scvi.model.SCVI(adata)


# In[14]:



# Train the model.

# In[15]:


model.train()


# Save the model

# In[16]:


model.save(dir_checkpoint + "scvi_model_test/")


# In[10]:


#model = scvi.model.SCVI.load("scvi_models/", adata=adata, use_gpu=True)


# ### scVI Output

# In[17]:


latent = model.get_latent_representation()


# We can save the output of the model in the anndata object:

# In[18]:


adata.obsm["X_scVI"] = latent



adata.write_h5ad(dir_checkpoint + "anndata_object/")
