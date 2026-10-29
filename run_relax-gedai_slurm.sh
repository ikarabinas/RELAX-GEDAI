#!/bin/bash
#SBATCH --job-name=relax-gedai
#SBATCH --array=1-26
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/relax-gedai_%A_%a.log
#SBATCH --error=logs/relax-gedai_%A_%a.error

# Run RELAX-GEDAI processing in parallel by participant ID
# Specify the number of participant IDs to process in the above sbatch array (e.g. "1-10" for 10 participants, "%10" max processed at a time)
# Run this script with "sbatch run_relax-gedai_slurm.sh"
# The ICA method (and with it the output folder) comes from RELAX_ICA_METHOD, picard when unset:
#   sbatch --export=ALL,RELAX_ICA_METHOD=amica --job-name=relax-amica run_relax-gedai_slurm.sh

module load matlab/R2023a

# Export task ID
export SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID

# Run MATLAB and let it read the env var via getenv()
matlab -nodisplay -nosplash -nodesktop -r "run('RELAX_SET_PARAMETERS_AND_RUN_slurm.m'); exit"
