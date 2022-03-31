#!/usr/bin/env python
# coding: utf-8

# # Comparing embeddings

# In the following I will compare the latent spaces obtained when using either scRNA-seq, scATAC-seq or both modalities to create the embeddings. For each modality I tried 3 different  of latent space dimensions, namely, 10, 15, 20.

# In[1]:


import scanpy as sc
import scvi
import numpy as np


# In[2]:


from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_mutual_info_score


# In[3]:


import pandas as pd


# In[4]:


# create a dictionary to save ARI scores for different embeddings
metrics = {}


# In[5]:


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


# # scRNA-seq

# ## PCA

# In[6]:


adata = scvi.data.read_h5ad("gpu_trained_20_dim/anndata_object")


# In[7]:


# compute pca
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.1)


# In[8]:


# leiden clustering
sc.tl.leiden(adata, key_added="leiden_clusters", resolution=1.2)
print(f"Using Leiden clustering we obtain: {len(adata.obs.leiden_clusters.unique())} clusters")


# In[9]:


# kmeans clustering
emb = adata.obsm["X_pca"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)
adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[10]:


sc.pl.umap(adata, color=["leiden_clusters", "kmeans", "celltype"], title=["leiden - pca", "kmeans", "celltype"],palette = colPalette_celltypes, size = 15)


# In[11]:


metrics["pca"] =  (adjusted_rand_score(adata.obs.celltype, adata.obs.leiden_clusters),adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans),
               adjusted_mutual_info_score(adata.obs.celltype, adata.obs.leiden_clusters), 
               adjusted_mutual_info_score(adata.obs.celltype, adata.obs.kmeans))


# ## scVI

# In[12]:


paths = ["gpu_trained_10_dim/", "gpu_trained_15_dim/", "gpu_trained_20_dim/"]


# In[13]:


for p in paths:
    adata = scvi.data.read_h5ad(p + "anndata_object")
    model = scvi.model.SCVI.load(p, adata=adata, use_gpu=False)
    # use scVI latent space for UMAP generation
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.3)
    
    # leiden clustering
    sc.tl.leiden(adata, key_added="scVI_clusters", resolution=1.2)
    print(f"Using Leiden clustering we obtain: {len(adata.obs.scVI_clusters.unique())} clusters")
    print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters)}")
    
    # kmeans clustering
    emb = adata.obsm["X_scVI"]
    kmeans = KMeans(n_clusters=37,
                    init="random",
                    n_init=200,
                    random_state=0).fit(emb)
    adata.obs["kmeans"] = kmeans.labels_.astype(str)
    
    sc.pl.umap(adata, color=["scVI_clusters", "kmeans", "celltype"], title = [f"scVI_cluster - {p}", "kmeans", "celltype"], palette = colPalette_celltypes, size = 15)

    metrics[p] =  (adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters),adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans),
                   adjusted_mutual_info_score(adata.obs.celltype, adata.obs.scVI_clusters), 
                   adjusted_mutual_info_score(adata.obs.celltype, adata.obs.kmeans))


# In[ ]:





# In[ ]:





# In[ ]:





# # scATAC-seq

# In[14]:


paths = ["atac_gpu_trained_10_dim/", "atac_gpu_trained_15_dim/", "atac_gpu_trained_20_dim/"]


# In[ ]:





# In[16]:


for p in paths:
    adata = scvi.data.read_h5ad(p + "anndata_object_peakvi")
    model = scvi.model.PEAKVI.load(p, adata=adata, use_gpu=False)
    # use scVI latent space for UMAP generation
    sc.pp.neighbors(adata, use_rep="X_PeakVI")
    sc.tl.umap(adata, min_dist=0.3)
    
    # leiden clustering
    sc.tl.leiden(adata, key_added="peakVI_clusters", resolution=1.2)
    print(f"Using Leiden clustering we obtain: {len(adata.obs.peakVI_clusters.unique())} clusters")
    
    # kmeans clustering
    emb = adata.obsm["X_PeakVI"]
    kmeans = KMeans(n_clusters=37,
                    init="random",
                    n_init=200,
                    random_state=0).fit(emb)
    adata.obs["kmeans"] = kmeans.labels_.astype(str)
    
    # plot umaps
    sc.pl.umap(adata, color=["peakVI_clusters", "kmeans", "celltype.mapped_mnn"], title = [f"peakVI_clusters - {p}", "kmeans", "celltype"], palette = colPalette_celltypes, size = 15)

    #print(f"The adjusted Rand index when using kmeans clustering is: {adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.peakVI_clusters)}")
    metrics[p] =  (adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.peakVI_clusters),adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.kmeans),
                   adjusted_mutual_info_score(adata.obs["celltype.mapped_mnn"], adata.obs.peakVI_clusters), 
                   adjusted_mutual_info_score(adata.obs["celltype.mapped_mnn"], adata.obs.kmeans))


# In[ ]:





# In[17]:


df = pd.DataFrame.from_dict(metrics)


# In[18]:


df = df.rename(index={0:"ARI-Leiden", 1:"ARI-kmeans", 2:"Mutual_information-Leiden", 3:"Mutual_information-kmeans"})


# In[19]:


df["metric"] = df.index
df.index = [0,1,2,3]


# In[ ]:


#df = df.rename(columns = {"gpu_trained_10_dim/":"model_rna_10dim", "gpu_trained_15_dim/":"model_rna_15dim", "gpu_trained_20_dim/":"model_rna_20dim", 
         #            "pca":"model_pca",
          #          "atac_gpu_trained_10_dim/":"model_atac_10dim", "atac_gpu_trained_15_dim/":"model_atac_15dim", "atac_gpu_trained_20_dim/":"model_atac_20dim"})


# In[ ]:


#df.columns = 'model_atac_10dim', 'model_atac_15dim', 'model_atac_20dim', 'model_pca','model_rna_10dim', 'model_rna_15dim', 'model_rna_20dim', 'metric'


# In[20]:


df


# In[21]:


df = df.melt("metric", var_name = "model", value_name = "val")


# In[22]:


import seaborn as sns
import matplotlib.pyplot as plt


# In[23]:


sns.barplot(y="model", x="val", hue="metric", data=df)
plt.legend(loc="lower left")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### 10 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("gpu_trained_10_dim/anndata_object")


# In[ ]:


model = scvi.model.SCVI.load("gpu_trained_10_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="scVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.scVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="scVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters)}")


# In[ ]:


metrics["rna_10_leiden"] =  adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_scVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["rna_10_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# ### 15 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("gpu_trained_15_dim/anndata_object")


# In[ ]:


model = scvi.model.SCVI.load("gpu_trained_15_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="scVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.scVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="scVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters)}")


# In[ ]:


metrics["rna_15_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_scVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["rna_15_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# ### 20 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("gpu_trained_20_dim/anndata_object")


# In[ ]:


model = scvi.model.SCVI.load("gpu_trained_20_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="scVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.scVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="scVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.scVI_clusters)}")


# In[ ]:


metrics["rna_20_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_scVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["rna_20_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# ## scATAC-seq

# ### 10 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("atac_gpu_trained_10_dim/anndata_object_peakvi")


# In[ ]:


pvi = scvi.model.PEAKVI.load("atac_gpu_trained_10_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_PeakVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="PeakVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.PeakVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="PeakVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)}")


# In[ ]:


metrics["atac_10_leiden"] = adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_PeakVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["atac_10_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# ### 15 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("atac_gpu_trained_15_dim/anndata_object_peakvi")


# In[ ]:


pvi = scvi.model.PEAKVI.load("atac_gpu_trained_15_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_PeakVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="PeakVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.PeakVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="PeakVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)}")


# In[ ]:


metrics["atac_15_leiden"] = adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_PeakVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["atac_15_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)


# ### 20 dimensional embedding

# #### Leiden Clustering

# In[ ]:


# load anndata with latent space
adata = scvi.data.read_h5ad("atac_gpu_trained_20_dim/anndata_object_peakvi")


# In[ ]:


pvi = scvi.model.PEAKVI.load("atac_gpu_trained_20_dim/", adata=adata, use_gpu=False)


# In[ ]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_PeakVI")
sc.tl.umap(adata, min_dist=0.3)


# In[ ]:


sc.tl.leiden(adata, key_added="PeakVI_clusters", resolution=1.2)


# In[ ]:


print(f"Using Leiden clustering we obtain: {len(adata.obs.PeakVI_clusters.unique())} clusters")


# In[ ]:


sc.set_figure_params(figsize=(8,8))
sc.pl.umap(adata, color="PeakVI_clusters", palette = colPalette_celltypes, size = 15)


# In[ ]:


print(f"The adjusted Rand index when using Leiden clustering is: {adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)}")


# In[ ]:


metrics["atac_20_leiden"] = adjusted_rand_score(adata.obs["celltype.mapped_mnn"], adata.obs.PeakVI_clusters)


# #### K-means clustering

# In[ ]:


emb = adata.obsm["X_PeakVI"]
kmeans = KMeans(n_clusters=37,
                init="random",
                n_init=200,
                random_state=0).fit(emb)


# In[ ]:


adata.obs["kmeans"] = kmeans.labels_.astype(str)


# In[ ]:


print(f"The adjusted Rand index when using k-means clustering is: {adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)}")


# In[ ]:


metrics["atac_20_kmeans"] = adjusted_rand_score(adata.obs.celltype, adata.obs.kmeans)

