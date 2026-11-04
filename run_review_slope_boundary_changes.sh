#!/bin/bash
#SBATCH --job-name=slope-review
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/review_slope_boundary_changes_%j.log
# Spectra and topographies of the components whose muscle call changes with boundary-aware slopes.
#   sbatch run_review_slope_boundary_changes.sh [cleaned file name]
module load matlab/R2023a
cd /home/nga4004/gitRepos/RELAX-GEDAI/.claude/worktrees/gentler-ICA
CLEANED=${1:-m544_dlpfc_day5_stitched_RELAX_GEDAI.set}
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; addpath(pwd); \
    root = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26_2exp/rest_treat_joined/'; \
    review_slope_boundary_changes([root 'RELAXProcessed_gentlerICA_picard'], '$CLEANED', [root 'slope_boundary_changes_picard/' extractBefore('$CLEANED', '_RELAX')])"
