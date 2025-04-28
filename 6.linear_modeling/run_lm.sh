#!/bin/bash

conda activate timelapse_analysis_env

jupyter nbconvert --to script --output-dir=scripts/ notebooks/*.ipynb

cd scripts || exit

python 0.linear_modeling.py

conda deactivate ; conda activate R_timelapse_env

Rscript 1.plot_linear_model_coefficients.r

cd ../ || exit

echo "Linear modeling and plotting completed."
