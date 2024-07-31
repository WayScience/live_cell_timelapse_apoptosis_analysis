#!/bin/bash

# this script is used to run the EDA process

# activate the conda environment
mamba activate timelapse_env

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts notebooks/*.ipynb

cd scripts

run the EDA script
python 0.generate_umap_embeddings.py --data_mode "CP"
python 0.generate_umap_embeddings.py --data_mode "scDINO"
python 0.generate_umap_embeddings.py --data_mode "combined"

# deactivate the conda environment
mamba deactivate

mamba activate R_timelapse_env

Rscript 1.visualize_umaps.r --data_mode "CP"
Rscript 1.visualize_umaps.r --data_mode "scDINO"
Rscript 1.visualize_umaps.r --data_mode "combined"

mamba deactivate

cd ../

# end of script
echo "EDA process completed"
