#!/bin/bash

conda init bash
conda activate timelapse_analaysis_env

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

cd scripts || exit

python 0.terminal_timepoint_anova.py
python 1.run_mAP_on_image_level_profiles.py

conda deactivate

conda activate R_timelapse_env

Rscript 2.plot_terminal_anova.r
Rscript 3.plot_mAP.r

conda deactivate
cd ../ || exit

echo "Completed the ground truth ananlysis"

