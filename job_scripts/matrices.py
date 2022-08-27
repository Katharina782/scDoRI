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


# create paths to read and write
dir_data = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/"
dir_checkpoint = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/"


#archr_gene_expr = scvi.data.read_h5ad(dir_data + "ArchR_gene_expr.h5ad")
archr_gene_scores_p2g = pd.read_csv(dir_data + "archr_gene_scores_p2g_table.csv", sep = " ", header=0)
archr_gene_scores = pd.read_csv(dir_data + "archr_gene_scores_table.csv", sep = " ", header=0)

# save gene names
#gene_list = archr_gene_scores_p2g.index.tolist()
#file_name = "gene_names.pkl"
#open_file = open(file_name, "wb")
#pickle.dump(gene_list, open_file)
#open_file.close()

# subset gene scores
#scores = archr_gene_scores[archr_gene_scores.index.isin(archr_gene_scores_p2g.var["name"])]
scores = archr_gene_scores.loc[archr_gene_scores_p2g.index.tolist(), :]


# convert to numpy array
scores_mat = scores.to_numpy()
#p2g_scores = archr_gene_scores_p2g.to_numpy()

# save as h5 file
h5f = h5py.File(dir_data + 'gene_scores.h5', 'w')
h5f.create_dataset('matrix', data=scores_mat)
h5f.close()

#h5f = h5py.File(dir_data + 'gene_scores_p2g.h5', 'w')
#h5f.create_dataset('matrix', data=p2g_scores)
#h5f.close()