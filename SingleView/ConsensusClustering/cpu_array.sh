#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3
#SBATCH --output=job.%A.%a.out

view=('clinical' 'contextual' 'physio' 'proteome')
p1=${view[`expr $SLURM_ARRAY_TASK_ID % ${#view[@]}`]}

python SingleView_DBSCAN_search.py --view $p1