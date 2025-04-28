#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --output=sc_sampling_parent-%j.out

module load miniforge
conda init bash
conda activate timelapse_map_env

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

# read percentage of parent cells from file
mapfile -t percentage < ./combinations/percentage.txt
mapfile -t seeds < ./combinations/seeds.txt

for p in "${percentage[@]}"; do
    for s in "${seeds[@]}"; do
            echo "Parent cell sampling: $p, seed: $s"
            number_of_jobs=$(squeue -u $USER | wc -l)
            while [ $number_of_jobs -gt 990 ]; do
                sleep 1s
                number_of_jobs=$(squeue -u $USER | wc -l)
            done
            sbatch child_sc_sampling_HPC.sh $p $s
    done
done

conda deactivate

echo "Parent cell sampling jobs submitted"
