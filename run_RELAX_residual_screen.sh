#!/bin/bash
#SBATCH --job-name=relax-residual-screen
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/residual_screen_%j.log
#SBATCH --error=logs/residual_screen_%j.error

# Residual-noise screen (RELAX_residual_screen.m): which cleaned files are clean enough for source localization.
#   mkdir -p logs && sbatch run_RELAX_residual_screen.sh [stage] [files]
#     stage - score | sheets | all (default: all)
#              sheets needs an earlier score run
#     files - optional comma-separated base names; sheets draws only these
#             e.g. m144_dlpfc_day1_stitched,m208_dlpfc_day1_stitched

STAGE="${1:-all}"
FILES="${2:-}"

mkdir -p logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath(pwd); RELAX_add_vendored_noisetools(); \
files = strsplit('${FILES}', ','); files = files(~cellfun(@isempty, files)); \
RELAX_residual_screen('${STAGE}', files)"
