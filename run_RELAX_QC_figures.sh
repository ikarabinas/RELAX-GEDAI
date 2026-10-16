#!/bin/bash
#SBATCH --job-name=relax-qc-figures
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/qc_figures_%j.log
#SBATCH --error=logs/qc_figures_%j.error

# Artifact-removal QC figure set (RELAX_QC_figures.m).
#   mkdir -p logs && sbatch run_RELAX_QC_figures.sh [stage] [files]
#     stage - extract | draw | all (default: all)
#     files - optional comma-separated base names to trial on a few files,
#             e.g. m144_dlpfc_day1_stitched,m208_dlpfc_day1_stitched

STAGE="${1:-all}"
FILES="${2:-}"

mkdir -p logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath(pwd); \
files = strsplit('${FILES}', ','); files = files(~cellfun(@isempty, files)); \
RELAX_QC_figures('${STAGE}', files)"
