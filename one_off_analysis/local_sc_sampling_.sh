#!/bin/bash

conda activate timelapse_map_env

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

# read percentage of parent cells from file
mapfile -t number_of_cells_per_well < ./combinations/number_of_cells.txt

cd scripts/ || exit
# get total length of percentage array
total=${#number_of_cells_per_well[@]}

iterator=0
for number_of_cells in "${number_of_cells_per_well[@]}"; do

    python 2.run_map_on_subsampled_cells.py --number_of_cells "$number_of_cells"
    python 2.run_map_on_subsampled_cells.py --number_of_cells "$number_of_cells" --shuffle
    iterator=$((iterator + 1))
    progress=$((iterator * 100 / $total))
    echo "Progress: $progress%"

done

cd .. || exit

conda deactivate

echo "Parent cell sampling jobs submitted"
