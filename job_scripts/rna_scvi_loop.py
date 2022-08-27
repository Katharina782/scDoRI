import scanpy as sc


# In[2]:


import scvi


# In[3]:


import numpy as np


# # RNA

# In[4]:

dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/scvi_training/"
anndata_object = "rna_anndata_from_seurat"
adata = scvi.data.read_h5ad(dir_data + anndata_object)






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
    n_top_genes=5000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="sample"
)


# In[12]:

for n in [10, 15, 20]:

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        categorical_covariate_keys=["sample"]
    )


    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",

        dropout_rate=0.2,
        n_layers=2,
        n_latent=n,
        gene_likelihood='nb'
    )

    model = scvi.model.SCVI(
        adata,
        **arches_params
    )



    # Train the model.


    model.train()


    # Save the model


    model.save(dir_checkpoint + f"{n}scvi_model/")


    latent = model.get_latent_representation()


    # We can save the output of the model in the anndata object:


    adata.obsm[f"X_scVI_{n}"] = latent

adata.write_h5ad(dir_checkpoint + "anndata_object_dims/")

    
for i in [1000, 5000, 10000]:
    
    adata = scvi.data.read_h5ad(dir_data + anndata_object)

    # Remove genes which are exprssed at very low levels or other outliers. Scvi tools use the raw count matrix. We can save the raw count matrix in an extra layer and then normalize the count matrix for other tasks. 
    sc.pp.filter_genes(adata, min_counts=3)


    adata.layers["counts"] = adata.X.copy() # preserve counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata # freeze the state in `.raw`


    # Reduce the number of features which will be used as input for the model. You can use 1,000 to 10,000 most highly variable genes. Lets try with 10,000 first. Here we can already add a batch_key which means that the highly variable genes will be selected for each "batch" separately. We will not select all highly variable genes from just one "batch". In our case the only batch effect we would expect is between cells of different replicates. 

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=i,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="sample"
    )

    
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        categorical_covariate_keys=["sample"]
    )


    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",

        dropout_rate=0.2,
        n_layers=2,
        n_latent=15,
        gene_likelihood='nb'
    )

    model = scvi.model.SCVI(
        adata,
        **arches_params
    )



    # Train the model.


    model.train()


    # Save the model


    model.save(dir_checkpoint + f"{i}_hvg_15dim_scvi_model/")


    latent = model.get_latent_representation()


    # We can save the output of the model in the anndata object:


    adata.obsm[f"X_scVI_{i}_hvg"] = latent

    adata.write_h5ad(dir_checkpoint + f"anndata_object_{i}hvg/")

    