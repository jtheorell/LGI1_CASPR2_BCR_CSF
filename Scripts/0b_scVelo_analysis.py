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

loompy.combine(["onefilepercell_JR1166_1_1_B02_and_others_N8PVQ.loom","onefilepercell_JR1227_1_1_A05_and_others_LBT7E.loom","onefilepercell_JR1284_1_1_A12_and_others_O2FQT.loom"], "all_donors.loom")

scv.set_figure_params()

adata = scv.read("all_donors.loom")

scv.pl.proportions(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

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
