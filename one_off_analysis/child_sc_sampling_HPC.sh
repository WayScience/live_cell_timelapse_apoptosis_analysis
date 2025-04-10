#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --output=sc_sampling_child-%j.out

module load miniforge
conda init bash
conda activate timelapse_map_env

percent=$1
seed=$2

cd scripts/ || exit


python 2.run_map_on_percentages_of_cells.py --percentage $percent --seed $seed --shuffle
python 2.run_map_on_percentages_of_cells.py --percentage $percent --seed $seed

cd .. || exit

conda deactivate

echo "Child cell sampling job completed"
