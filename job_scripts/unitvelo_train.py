import unitvelo as utv
import scvelo as scv
scv.settings.verbosity = 0
import scanpy as sc
#import tensorflow as tf



dir_data = "/dkfz/cluster/gpu/data/OE0533/k552k/"
dir_checkpoint = "/dkfz/cluster/gpu/checkpoints/OE0533/k552k/"


velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.GPU = 0



label = 'celltype.mapped'


adata = utv.run_model(dir_data + "mesoderm.h5ad", label, config_file=velo_config)



sc.write(dir_checkpoint + "unitvelo_mesoderm", adata)
