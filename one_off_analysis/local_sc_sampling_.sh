#!/bin/bash

conda activate timelapse_map_env

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

# read percentage of parent cells from file
mapfile -t percentage < ./combinations/percentage.txt
mapfile -t seeds < ./combinations/seeds.txt

cd scripts/ || exit
# get total length of percentage array
total=${#percentage[@]}
# get total length of seeds array
total_seed=${#seeds[@]}

# get the total length of the array
total_length=$((total * total_seed))
iterator=0
for percent in "${percentage[@]}"; do
    for seed in "${seeds[@]}"; do
        python 2.run_map_on_percentages_of_cells.py --percentage "$percent" --seed "$seed"
        python 2.run_map_on_percentages_of_cells.py --percentage "$percent" --seed "$seed" --shuffle
        iterator=$((iterator + 1))
        progress=$((iterator * 100 / total_length))
        echo "Progress: $progress%"
    done
done

cd .. || exit

conda deactivate

echo "Parent cell sampling jobs submitted"
