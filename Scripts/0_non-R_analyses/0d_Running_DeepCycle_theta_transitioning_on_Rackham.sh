ssh -AX jakob@rackham.uppmax.uu.se
#password is simply the standard one without the changing code
mkdir 230216_deepcycle_analysis
exit

cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/

rsync -ah DeepCycle jakob@rackham.uppmax.uu.se:/home/jakob/230216_deepcycle_analysis

ssh -AX jakob@rackham.uppmax.uu.se

cd 230216_deepcycle_analysis/DeepCycle

interactive -A snic2022-22-1086

module load conda
#conda env create -f DeepCycle_env.yml
conda activate DeepCycle

pip install -upgrade anndata

cd 230216_deepcycle_analysis/DeepCycle

python ./scripts/estimate_cell_cycle_transitions.py \
  --input_adata ./adata_DeepCycle.h5ad \
  --gene_phase_dict ./theta_annotation/gene_phase_dict.json \
  --output_npy_transitions ./theta_transitions_fibro.npy \
  --output_svg_plot ./cell_cycle_phase_detection_fibro.svg

sbatch transition_definition.sh

jobinfo -u jakob

#So this does not work, as we have too few genes associated to the cell cycle expressed. No CyclinE genes are expressed, e.g. 
