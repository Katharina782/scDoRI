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


# create paths to read and write
dir_data = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/"
dir_checkpoint = "/omics/groups/OE0533/internal/katharina/scDoRI/gastrulation_data/jupyter_notebooks/"


#archr_gene_expr = scvi.data.read_h5ad(dir_data + "ArchR_gene_expr.h5ad")
archr_gene_scores_p2g = pd.read_csv(dir_data + "archr_gene_scores_p2g_table.csv", sep = " ", header=0)
archr_gene_scores = pd.read_csv(dir_data + "archr_gene_scores_table.csv", sep = " ", header=0)

# save gene names
gene_list = archr_gene_scores_p2g.index.tolist()
file_name = "gene_names.pkl"
open_file = open(file_name, "wb")
pickle.dump(gene_list, open_file)
open_file.close()

# subset gene scores
#scores = archr_gene_scores[archr_gene_scores.index.isin(archr_gene_scores_p2g.var["name"])]
#scores = archr_gene_scores.loc[archr_gene_scores_p2g.index.tolist(), :]

