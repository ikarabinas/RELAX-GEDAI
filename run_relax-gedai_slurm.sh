#!/bin/bash
#SBATCH --job-name=relax-gedai
#SBATCH --array=1-33%17
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/relax-gedai_%A_%a.log
#SBATCH --error=logs/relax-gedai_%A_%a.error

# Run RELAX-GEDAI processing in parallel by participant ID
# Specify the number of participant IDs to process in the above sbatch array (e.g. "1-10" for 10 participants, "%10" max processed at a time)
# Run this script with "sbatch run_relax-gedai_slurm.sh"

module load matlab/R2023a

# Export task ID
export SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID

# Run MATLAB
matlab -nodisplay -nosplash -nodesktop -r "run('RELAX_SET_PARAMETERS_AND_RUN_slurm.m'); exit"
