#!/bin/bash
#SBATCH --job-name=muscle-compare
#SBATCH --partition=scu-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=0-06:00:00
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.error
# Runs compare_muscle_arms on a compute node after both run_arm.sh arms finish.
# Submit from this directory: sbatch run_compare.sh
set -eo pipefail
here=$(pwd)
worktree=$(cd .. && pwd)
module load matlab/R2023a
matlab -nodisplay -nosplash -nodesktop -batch "addpath('/home/nga4004/tools/eeglab2026.0.0/'); eeglab nogui; addpath('/home/nga4004/tools/eeglab2026.0.0/plugins/Fieldtrip-lite250523/'); addpath('$worktree/'); addpath('$here'); RELAX_add_vendored_noisetools(); compare_muscle_arms"
