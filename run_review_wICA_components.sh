#!/bin/bash
#SBATCH --job-name=review-wica
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/review_wICA_components_%A_%a.log
# Per-component review of the gentler-ICA PICARD and AMICA runs (review_wICA_components.m).
#   sbatch --array=1-89 run_review_wICA_components.sh        one task per cleaned file
#   sbatch --dependency=afterany:<array job> --export=ALL,REVIEW_MODE=combine run_review_wICA_components.sh
module load matlab/R2023a
cd /home/nga4004/gitRepos/RELAX-GEDAI/.claude/worktrees/gentler-ICA
MODE=${REVIEW_MODE:-extract}
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; addpath(pwd); \
    root = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26_2exp/rest_treat_joined/'; \
    runDirs = {[root 'RELAXProcessed_gentlerICA_picard'], [root 'RELAXProcessed_gentlerICA_amica']}; \
    outDir = [root 'wICA_components_picard_vs_amica']; \
    if strcmp('$MODE', 'combine'); review_wICA_components('combine', runDirs, {'picard', 'amica'}, outDir); \
    else; review_wICA_components('extract', runDirs, {'picard', 'amica'}, outDir, ${SLURM_ARRAY_TASK_ID:-1}); end"
