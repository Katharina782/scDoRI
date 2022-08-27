import scanpy as sc
import scvi
import numpy as np
import anndata as ad
from anndata import AnnData
from scipy import sparse
import pandas as pd

# paths to files and save
dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"

# read rna and atac anndata objects
rna_adata = scvi.data.read_h5ad(dir_data + "anndata_rna.h5ad")
atac_adata = scvi.data.read_h5ad(dir_data + "anndata_atac_peak_matrix.h5ad")

#adata.var_names_make_unique()
rna_adata.var_names_make_unique()
atac_adata.var_names_make_unique()

# this old object is used for filtering, because we will use those gene annotations
old = scvi.data.read_h5ad(dir_data + "old_anndata_rna.h5ad")

rna_adata = rna_adata[:,rna_adata.var.index.isin(old.var.index)]
annotations = pd.merge(old.var, rna_adata.var, on="gene", how = "inner")
annotations.index = old.var.index
annotations = annotations.drop(["gene", "Strand"], axis=1)
annotations = annotations.rename(columns={"Accession":"ID", "End":"end", "Start":"start", "Chromosome":"chr"})
rna_adata.var = annotations

# combine counts and accessibility matrix
matrix = sparse.hstack((rna_adata.X, atac_adata.X))

# create dataframe for rna
rna_df = rna_adata.var
rna_df["modality"] = ["Gene_expression" for i in range(len(rna_df))]

# create dataframe for atac
atac_df = atac_adata.var
atac_df = atac_df.rename(columns={"idx":"ID"})
atac_df["modality"]=["Peaks" for i in range(len(atac_df))]
atac_df = atac_df.drop(["score"], axis = 1)

# combine the two dataframes
df = rna_df.append(atac_df)

# create a combined anndata object
multiome = AnnData(
X = matrix,
obs = atac_adata.obs,
var = df,
uns = rna_adata.uns)

# create a multiome object
adata_mvi = scvi.data.organize_multiome_anndatas(multiome)

# sort the object var
adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()


# filter features which appear in < 1%
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))


# setup anndata for the model
scvi.model.MULTIVI.setup_anndata(adata_mvi, categorical_covariate_keys=['Sample'])#, batch_key='modality')


mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var['modality']=='Gene_expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
    dropout_rate=0.2,
    n_layers_encoder=2,
    n_layers_decoder=2,
    n_latent=15,
    latent_distribution="normal",
)


# train the model
mvi.train()




# save model
mvi.save(dir_checkpoint + "multivi_model/")

# extract latent embedding
latent = mvi.get_latent_representation()

# add latent embedding to anndata object
adata_mvi.obsm["X_multivi"] = latent


adata_mvi.write_h5ad(dir_checkpoint + "anndata_object/")
