import scanpy as sc
import scvi
import numpy as np
from scipy.sparse import csr_matrix


dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/peakvi_training/"
anndata_object = "peak_anndata_from_archr"
adata = scvi.data.read_h5ad(dir_data + anndata_object)


min_cells = int(adata.shape[0] * 0.05)

# filter the regions using scanpy
sc.pp.filter_genes(adata, min_cells=min_cells)


# convert sparse matrix to csr matrix for faster training
copy_mat = adata.X.copy()
csr_matrix = csr_matrix(copy_mat)
adata.layers["csr"] = csr_matrix



for n in [10, 15, 20]:

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
        n_latent=n,
        latent_distribution="normal"
    )



    pvi = scvi.model.PEAKVI(
        adata,
        **arches_params
    )


    # train the model
    pvi.train()


    # save the model
    pvi.save(dir_checkpoint + f"peakvi_model_{n}_dim/")

    # get latent space
    latent = pvi.get_latent_representation()
    adata.obsm[f"X_PeakVI_{n}_dim"] = latent


adata.write_h5ad(dir_checkpoint + "anndata_object_peakvi_dims/")

