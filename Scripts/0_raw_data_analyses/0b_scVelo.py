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

os.chdir("/Users/jakob.theorell/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/Velocity/velocyto_out")
#We also verify that we have ended up in the right place
# verify the path using getcwd()
cwd = os.getcwd()
 
# print the current directory
print("Current working directory is:", cwd)

loompy.combine(["onefilepercell_JR1166_1_1_B02_and_others_N8PVQ.loom","onefilepercell_JR1227_1_1_A05_and_others_LBT7E.loom","onefilepercell_JR1284_1_1_A12_and_others_O2FQT.loom"], "all_donors.loom")

scv.set_figure_params()

adata = scv.read("all_donors.loom")

scv.pl.proportions(adata)

#I have increased this to 4000 to make sure that we have some cyclin genes included for the cell cycle step. 
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=4000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#Here, we will export the adata file, for later use in DeepCycle. 
adata.write("../../DeepCycle/adata_post_moments.h5ad")

scv.tl.velocity(adata, mode='deterministic')

scv.tl.velocity_graph(adata)

scv.tl.umap(adata)
scv.tl.louvain(adata, flavor = 'igraph')
scv.pl.velocity_embedding_stream(adata, basis='umap')

scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

#And can that last thing be exported, i.e. the peudotime vector? Would be extremely useful. 
adata.write_csvs('All_cells')
adata.write("../../DeepCycle/adata_complete.h5ad")

quit()

