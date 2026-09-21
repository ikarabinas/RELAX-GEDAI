#!/bin/bash
#SBATCH --job-name=relax-refpop-survey
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/refpop_survey_%j.log
#SBATCH --error=logs/refpop_survey_%j.error

# Real-data survey of reference-electrode pops on five files (survey_reference_pops.m).
# Submit with: mkdir -p logs && sbatch run_survey_reference_pops.sh
# Writes reference_pop_survey.mat in the repo root and PNGs under
# docs/superpowers/reports/reference_pop_survey/.

mkdir -p logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); \
addpath(pwd); \
survey_reference_pops('reference_pop_survey.mat');"
