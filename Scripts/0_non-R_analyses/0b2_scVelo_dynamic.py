#This analysis is much dependent on 
#https://scvelo.readthedocs.io/getting_started/

#Installtion: 
#pip install -U scvelo
#pip install python-igraph louvain
#pip install pybind11 hnswlib

cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/Velocity/velocyto_out
python3
import scvelo as scv
import loompy
import os
import pandas as pd
import numpy as np
import scanpy as sc

os.chdir("/Users/jakob.theorell/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/Velocity/velocyto_out")
#We also verify that we have ended up in the right place
# verify the path using getcwd()
cwd = os.getcwd()
 
# print the current directory
print("Current working directory is:", cwd)

loompy.combine(["onefilepercell_JR1166_1_1_B02_and_others_N8PVQ.loom","onefilepercell_JR1227_1_1_A05_and_others_LBT7E.loom","onefilepercell_JR1284_1_1_A12_and_others_O2FQT.loom"], "all_donors.loom")

scv.set_figure_params()

adata = scv.read("all_donors.loom")

#Now, the cells that should be included has been somewhat reduced since previously. Therefore, we need to export a list of the cell names, so that the correct ones can be included. 
adata.write_csvs('Raw_cells')

#Now, we re-import the data that will tell us which cells to include
colData = pd.read_csv("Raw_cells/obs_plus_cell_type.csv")
adata.obs = colData
adata = adata[adata.obs.excluded == False]


#We are now depending heavily on this document: 
#https://scvelo.readthedocs.io/en/stable/DynamicalModeling/

scv.pl.proportions(adata)

#We are going with a high number of genes, as we have an unusually deep dataset and as we are interested in cell cycle dynamics, which are not detectable among the most highly variable genes. 
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=12000)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#Now, here comes a new thing: the addition of dynamics
scv.tl.recover_dynamics(adata, n_jobs = 8)
scv.tl.velocity(adata, mode='dynamical')

scv.tl.velocity_graph(adata)

scv.tl.umap(adata)
scv.tl.louvain(adata, flavor = 'igraph')

scv.tl.latent_time(adata)
scv.tl.score_genes_cell_cycle(adata)

scv.tl.velocity_confidence(adata)

scv.tl.velocity_pseudotime(adata)

#Now some plotting: 
#adata = scv.read("adata_complete.h5ad")

scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'singler', palette = ["#31688EFF", "#FDE725FF", "#440154FF", "#35B779FF"], alpha = 0.7, add_margin = 0.05, size = 300, save = "velocity_embedding.svg", figsize = (5,5))

scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', alpha = 0.7, figsize = (5,5), save = "latent_time.svg") #size = 80, 
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], figsize = (5,5), save = "cell_cycle.svg") #alpha = 0.9, 
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], figsize = (5.5,5), save = "confidence.svg")
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', figsize = (5,5), save = "pseudotime.svg")

adata.obs["donor"] = [x[17:21] for x in adata.obs.CellID]
adata.obs["antigen"] = [1 if item == '1166' else 1 if item == '1227' else 2 if item == "1284" else item for item in adata.obs.donor]
adata.obs["antigen"] = pd.Categorical(adata.obs["antigen"])
scv.pl.scatter(adata, color=adata.obs["antigen"], palette = ["#FAA51A", "#FF6633"], figsize = (5,5), save = "antigen.svg", alpha = 0.2)
#This seems to indicate a perfect separation between the groups in some places, but this is not really a problem, as the larger PC-containgin region as well as the B-cell region spans both groups. 

#And can that last thing be exported, i.e. the peudotime vector? Would be extremely useful. 
adata.write_csvs('All_cells')
adata.write("adata_complete.h5ad")
