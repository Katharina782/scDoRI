import scanpy as sc
import scvi
import numpy as np
import anndata as ad
from anndata import AnnData
from scipy import sparse


# create paths to read and write
dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"
adata = scvi.data.read_h5ad(dir_data + "anndata_rna.h5ad")

# read in files
rna_adata = scvi.data.read_h5ad(dir_data + "anndata_rna.h5ad")
atac_adata = scvi.data.read_h5ad(dir_data + "anndata_atac_peak_matrix.h5ad")

# create a joint matrix

matrix = sparse.hstack((rna_adata.X, atac_adata.X))

# create a dataframe for the var slot of the new anndata object

# rna
rna_df = rna_adata.var
rna_df = rna_df.rename(columns = {"gene":"ID"})
rna_df["modality"] = ["Gene_expression" for i in range(len(rna_df))]

# atac
atac_df = atac_adata.var
atac_df = atac_df.rename(columns={"chr":"ID"})
atac_df["modality"]=["Accessibility" for i in range(len(atac_df))]

# combine rna & atac
df = rna_df.append(atac_df)



# create anndata object
multiome = AnnData(
X = matrix,
obs = atac_adata.obs,
var = df,
uns = rna_adata.uns)


# prepare for setting up the model
#multiome_mvi = scvi.data.organize_multiome_anndatas(multiome)


# set
scvi.model.MULTIVI.setup_anndata(multiome_mvi,
                                 batch_key = "modality",
                                 categorical_covariate_keys=["sample"])
                                 
# create model
arches_params = dict(
    n_genes=(multiome_mvi.var["modality"]=="Gene_expression").sum(),
    n_regions=(multiome_mvi.var["modality"]=="Peaks").sum(),
    dropout_rate=0.2,
    n_layers_encoder=2,
    n_layers_decoder=2,
    n_latent=20,
    latent_distribution='normal'
)

mvi = scvi.model.MULTIVI(
    multiome_mvi,
    **arches_params
)

# train the model
mvi.train()

# save model
mvi.save(dir_checkpoint + "multivi_model/")

# extract latent embedding
latent = mvi.get_latent_representation()

# add latent embedding to anndata object
multiome.obsm["X_multivi"] = latent


multiome.write_h5ad(dir_checkpoint + "anndata_object/")
