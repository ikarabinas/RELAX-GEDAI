#!/bin/bash
#SBATCH --job-name=compare-wica
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/compare_wICA_runs_%j.log
# Paired comparison of the gentler-ICA PICARD baseline against the AMICA run (compare_wICA_runs.m).
module load matlab/R2023a
cd /home/nga4004/gitRepos/RELAX-GEDAI/.claude/worktrees/gentler-ICA
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; addpath(pwd); \
    root = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26_2exp/rest_treat_joined/'; \
    compare_wICA_runs({[root 'RELAXProcessed_gentlerICA_picard'], [root 'RELAXProcessed_gentlerICA_amica']}, \
        {'picard', 'amica'}, [root 'wICA_comparison_picard_vs_amica'])"
