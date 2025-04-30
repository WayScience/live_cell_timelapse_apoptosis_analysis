#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=00:5:00
#SBATCH --output=sc_sampling_child-%j.out

module load miniforge
conda init bash
conda activate timelapse_map_env

number_of_cells=$1

cd scripts/ || exit


python 2.run_map_on_subsampled_cells.py --number_of_cells "$number_of_cells"
python 2.run_map_on_subsampled_cells.py --number_of_cells "$number_of_cells" --shuffle

cd .. || exit

conda deactivate

echo "Child cell sampling job completed"
