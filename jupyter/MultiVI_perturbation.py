

import scanpy as sc
import scvi
import numpy as np
import anndata as ad
from anndata import AnnData
from scipy import sparse
import pandas as pd





# create paths to read and write
dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"

rna_adata =  sc.read_h5ad(data_dir + "rna_anndata_from_seurat")
atac_adata = sc.read_h5ad(data_dir + "peak_anndata_from_archr")




# make sure that the cells are in the same order
atac_adata = atac_adata[rna_adata.obs.index.tolist(), :]


atac_df = atac_adata.var
atac_df["modality"]=["Peaks" for i in range(len(atac_df))]





rna_df = rna_adata.var
#rna_df = rna_df.rename(columns = {"gene":"ID"})
rna_df["modality"] = ["Gene_expression" for i in range(len(rna_df))]



df = rna_df.append(atac_df)



matrix = sparse.hstack((rna_adata.X, atac_adata.X))



multiome = AnnData(
X = matrix,
obs = atac_adata.obs,
var = df)



adata_mvi = scvi.data.organize_multiome_anndatas(multiome)



adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()



scvi.model.MULTIVI.setup_anndata(adata_mvi)



mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var['modality']=='Gene_expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
)


mvi.train()


# save model
mvi.save(dir_checkpoint + "multivi_model/")

# extract latent embedding
latent = mvi.get_latent_representation()

# add latent embedding to anndata object
multiome.obsm["X_multivi"] = latent


multiome.write_h5ad(dir_checkpoint + "anndata_object/")