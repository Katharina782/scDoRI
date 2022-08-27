#!/usr/bin/env python
# coding: utf-8

# # Creating cell aggregates

# We want to create 500 cell aggreagates of 100 cells each. To create these groupings, we will randomly sample 500 cell from the latent space embedding. Using a nearest neighbor search we will then find the 100 nearest neighbors of these cells. Any group of cells which has an overlap >80% with any of the previously created groups will be removed. These cell aggregates can then be used for downstream correlation tasks.

# In the following I subset all matrices (gene expression, ArchR gene scores computed from scATAC-seq and gene scores computed from peak-to-gene links to contain an overlapping list of genes. The reason why each matrix has a different number of genes is that ArchR gene scores contain more genes, because based on chromatin accessibility we might 

# In[8]:


import scanpy as sc
import scvi
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import NearestNeighbors
from numba import jit, njit
import seaborn as sns
from scipy import sparse
import seaborn as sns
import pandas as pd
import h5py
import pickle


dir_data = "/omics/groups/OE0533/internal/katharina/scDori/gastrulation/jupyter_notebooks/"


# In[23]:


gene_names = pd.read_pickle(r"/omics/groups/OE0533/internal/katharina/scDori/gastrulation/jupyter_notebooks/gene_names.pkl")


# In[26]:


# read in the anndata object which contains the latent space embedding
adata = scvi.data.read_h5ad(dir_data + "gpu_trained_20_dim/anndata_object")

# read the anndata object containing ArchR gene expression matrix
archr_gene_expr = scvi.data.read_h5ad(dir_data + "ArchR_gene_expr.h5ad")
# subset anndata object to contain only overlapping cells
adata = adata[archr_gene_expr.obs.index, :]

# get latent space embedding
latent_embedding = adata.obsm["X_scVI"]


# Prepare the metadata for the cell aggregate function:

# In[28]:


metadata = adata.obs
# create column for cell names
metadata["cells"] = metadata.index
# create index for cells
metadata["idx"] = np.arange(len(metadata))


# https://stackoverflow.com/questions/20928136/input-and-output-numpy-arrays-to-h5py

# ### Sample cells & compute nearest neighbors
# 
# Since we have 45,991 cells in our dataset, we will sample 1000 cell aggregates of 50 cells each. 

# In[29]:


class nearest_neighbors:
    def __init__(self, latent_embedding):
    # attribute root of the class Tree will be an instance of the class Node
    # attriute self.root is an object of class Node
        self.metadata = metadata
        self.latent_embedding = latent_embedding
        
    def sampling_cells(self, n_aggregates):
        self.sample_cells = np.random.choice(self.latent_embedding.shape[0], n_aggregates, replace=False)
    
    def compute_NN(self, k):
        nbrs = NearestNeighbors(n_neighbors=k, algorithm="ball_tree").fit(self.latent_embedding)
        dist, ind = nbrs.kneighbors(self.latent_embedding)
        # subset the nearest neighbors to only contain the sampled cells
        self.distance = dist[self.sample_cells, ]
        self.index = ind[self.sample_cells, ]

        #print(f"Out of {nrow} cell aggregates we keep {keep.shape} cell aggregates which do not overlap with another aggregate by more than 80%.")
        # filter the matrix for cells to keep
        #self.idx_keep = self.index[keep, :]
    


# ### Check for overlapping cell aggregates

# In[30]:


# We want to check whether the 50 nearest neighbors of any given cell overlap with the 50 nearest neighbors of any other cell by more than 80%
@jit(nopython=True)
def check_overlap(index):
    nrow = index.shape[0]
    ncol = index.shape[1]
    # create an array to store whether a cell passed the overlap check
    # all entries are initially False
    considered_cells = np.zeros(nrow).astype(np.bool8)
    # loop over each cell and cosider to add it to the set of cell aggregates
    for i in range(nrow):
        check = True
        # loop over all previous aggregates 
        for comp in np.where(considered_cells)[0]:
            # get the number of cells which overlap between the neighborhood of the cell we would like to add and the neighborhoud of previous cell "comp"
            intersect = np.intersect1d(index[i, :], index[comp, :])
            # for each comparison between current cell i which we would like to add and previous cell which we are comparing to
            # compute the percentage of overlap
            if (len(intersect) / ncol) > 0.8: # if the intersection is above 0.8, we do not consider it
                check = False
                break
        if check:
            considered_cells[i] = True
    # get indices
    keep = np.arange(start=0, stop=nrow, step=1)[considered_cells]
    return index[keep, :]
    


# ### Only keep cells of same celltype

# In[31]:



# We only want to keep cells of the same celltype
def filter_celltypes(metadata, sample_cells, idx_keep):
    groups={}
    # check whether a cell aggregate contains cells of other celltypes and remove then
    for n, i in enumerate(sample_cells):
        # get celltype information of sampled cell
        celltype_test_cell = metadata.iloc[i]["celltype"]
        # get indices of cells which are in the neighborhood 
        neighbor_cells = idx_keep[idx_keep[:, 0] == i, :]
        # get cells which are of the same celltype, vector includes the sampled cell itself
        keep = np.array(metadata[(metadata.idx.isin(neighbor_cells.flatten())) & (metadata.celltype == celltype_test_cell)]["idx"])
        # keep only aggregates which contain at least 10 cells after removing non-matching celltypes
        if keep.shape[0] > 10:
            groups[i] = keep
        else:
            continue
    # add dictionary to self
    return groups


# #### Apply functions to our data

# In[32]:


agg_object = nearest_neighbors(latent_embedding) # initialize the cell aggreagte object


# In[33]:


agg_object.sampling_cells(n_aggregates=1000)


# In[34]:


agg_object.compute_NN(k=50) # compute nearest neighbors 


# In[35]:


idx_keep = check_overlap(agg_object.index)


# In[37]:


metadata = adata.obs
# create column for cell names
metadata["cells"] = metadata.index
# create index for cells
metadata["idx"] = np.arange(len(metadata))


# In[38]:


groups = filter_celltypes(metadata, agg_object.sample_cells,idx_keep)





# ## Aggregate expression/scores/accessibility 

# In[48]:


archr_gene_expr.var["index"] = archr_gene_expr.var.index


# In[50]:


new = archr_gene_expr.var[archr_gene_expr.var["name"].isin(gene_names)]


# In[66]:


index_list = pd.to_numeric(new["index"]).to_numpy()


# In[39]:


rna_mat = archr_gene_expr.X


# TODO: Implement aggregation with numba -> convert to csr matrx

# In[40]:


def create_aggregates(mat, groups):
    # initialize matrix to store average gene expression for each cell aggregate
    # the matrix has dimensions genes x cell aggregates
    rna_agg = np.zeros((mat.shape[1],len(groups)))
    # for each cell aggregate calculate the average gene expression for each gene
    for i, g in enumerate(groups):
        rna_agg[:, i] = mat[groups[g], :].mean(axis=0).flatten()
    return sparse.csr_matrix(rna_agg) # return sparse matrix genes x cell aggregates


# In[ ]:





# In[41]:


# gene expression
rna_agg = create_aggregates(rna_mat, groups)


# In[68]:

# subset to contain only the genes which are also included in scores matrices
rna_agg = rna_agg[index_list, :]


# In[69]:



# #### Gene Scores from ArchR

# In[9]:


h5f = h5py.File(dir_data + 'gene_scores.h5','r')
scores = h5f['matrix'][:]
h5f.close()


# In[ ]:


scores = sparse.csr_matrix(scores)


# In[ ]:


# gene scores
score_agg = create_aggregates(scores.T, groups)


# #### Gene Scores from p2g links

# In[11]:


h5f = h5py.File(dir_data + "gene_scores_p2g.h5", "r")
scores_p2g = h5f["matrix"][:]
h5f.close()


# In[13]:


scores_p2g = sparse.csr_matrix(scores_p2g)




scores_p2g_agg = create_aggregates(scores_p2g.T, groups)


# ## Correlations

# To compute row-wise correlations:

# In[ ]:


@jit(nopython=True)
def rowwise_correlation(A, B, N):
    correlations = []
    for i in range(A.row(1)):
        # center and scale row i of matrix A
        rowa = A.row(i)
        rowa -= np.mean(rowa)
        rowa /= np.std(rowa)
        
        # center and scale row i of matrix B
        rowb = B.row(i)
        rowb -= np.mean(rowb)
        rowb /= np.std(rowb)
        
        # compute correlation between row i of matrix A and B
        corr = np.mean(rowa*rowb)
        correlations.append(corr)
    return correlations


# Convert sparse matrices to csr matrices:

# In[ ]:


expr_agg = csr.CSR.from_scipy(rna_agg)


# In[ ]:


score_agg = csr.CSR.from_scipy(score_agg)


# In[ ]:


scores_p2g_agg = csr.CSR.from_scipy(scores_p2g_agg) 


# In[ ]:


# get number of genes
N = rna_agg.shape[0]


# In[ ]:


corr_expr_scores = rowwise_correlation(expr, scores, N)


# In[ ]:


corr_expr_p2g_scores = rowwise_correlation(rna_agg, p2g_score_agg, N)


corr_scores = rowwise_correlation(p2g_score_agg, scores, N)


# ### Save correlations

# In[ ]:


with open(dir_data + 'corr_expr_scores.pkl', 'wb') as f:
    pickle.dump(corr_expr_scores, f)


# In[ ]:


with open(dir_data + 'corr_expr_p2g_scores.pkl', 'wb') as f:
    pickle.dump(corr_expr_p2g_scores, f)


    
with open(dir_data + 'corr_scores.pkl', 'wb') as f:
    pickle.dump(corr_scores, f)





