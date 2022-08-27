import scanpy as sc
import scvi


dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"
adata = scvi.data.read_h5ad(dir_data + "old_anndata_rna.h5ad")

min_cells = int(adata.shape[0] * 0.05)

# filter the regions using scanpy
sc.pp.filter_genes(adata, min_cells=min_cells)

# convert sparse matrix to csr matrix for faster training
from scipy.sparse import csr_matrix
copy_mat = adata.X.copy()
csr_matrix = csr_matrix(copy_mat)
adata.layers["csr"] = csr_matrix

# prepare anndata for the model  
scvi.model.PEAKVI.setup_anndata(adata,  
                                layer = "csr",
                                categorical_covariate_keys=["sample"])

# define model parameters 
arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    
    dropout_rate=0.2,
    n_layers_encoder=2,
    n_layers_decoder=2,
    n_latent=20,
    latent_distribution="normal"
)

pvi = scvi.model.PEAKVI(
    adata,
    **arches_params
)


# train the model
pvi.train()


# save the model
pvi.save(dir_checkpoint + "peakvi_model/")

# get latent space
latent = pvi.get_latent_representation()
adata.obsm["X_PeakVI"] = latent


# save anndata object with the latent space
adata.write_h5ad(dir_checkpoint + "anndata_object_peakvi/")
