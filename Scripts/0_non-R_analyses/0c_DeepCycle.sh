
#Here, we will aim to run DeepCycle.

#The DeepCycle software can be found  here: 
#https://github.com/andreariba/DeepCycle
#It is from this publication: 
#https://www.nature.com/articles/s41467-022-30545-8#code-availability
#For this to work, we need to export the data we have generated

#After many attempts, it turned out that the most sensible strateghy for this was to run it on an intel-based mac, locally, 
#so this was done with the below script. 

#pip install tensorflow
#pip install tensorflow-probability 
#pip install tensorflow-datasets
cd ~/Desktop/DeepCycle
python3 DeepCycle.py \
  --input_adata adata_post_moments.h5ad \
  --gene_list go_annotation/GO_cell_cycle_annotation_human.txt \
  --base_gene NOTCH2 \
  --expression_threshold 0.5 \
  --output_adata adata_DeepCycle.h5ad

#After this, the file was re-imported into the current M1 environment. Then: 
cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/DeepCycle

python3

import scvelo as scv

adata = scv.read("adata_DeepCycle.h5ad")

adata.write_csvs('DeepCycle_csv_results')




