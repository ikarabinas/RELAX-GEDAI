#!/bin/bash
#SBATCH --job-name=slope-boundary
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/compare_slope_boundary_awareness_%j.log
# Whole-record versus per-segment muscle slopes on finished PICARD-run files.
#   sbatch run_compare_slope_boundary_awareness.sh
module load matlab/R2023a
cd /home/nga4004/gitRepos/RELAX-GEDAI/.claude/worktrees/gentler-ICA
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; addpath(pwd); \
    root = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26_2exp/rest_treat_joined/'; \
    compare_slope_boundary_awareness([root 'RELAXProcessed_gentlerICA_picard'], round(linspace(1, 89, 8)), [root 'slope_boundary_awareness_picard.csv'])"
