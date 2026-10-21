#!/bin/bash
#SBATCH --job-name=gedai-alpha
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=/home/nga4004/gitRepos/RELAX-GEDAI/logs/gedai_alpha_%A_%a.log
#SBATCH --error=/home/nga4004/gitRepos/RELAX-GEDAI/logs/gedai_alpha_%A_%a.error

# Broadband-GEDAI alpha figure pipeline (GEDAI_alpha_figures.m).
#   sbatch --array=1-N run_GEDAI_alpha_figures.sh extract
#     N is the session count printed by GEDAI_alpha_figures('list'); each array task
#     summarises one session and skips it if its cache file exists.
#   sbatch run_GEDAI_alpha_figures.sh STAGE
#     any other stage, as it becomes available.

STAGE="${1:?usage: sbatch [--array=1-N] run_GEDAI_alpha_figures.sh STAGE}"

mkdir -p /home/nga4004/gitRepos/RELAX-GEDAI/logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath(pwd); RELAX_add_vendored_noisetools(); GEDAI_alpha_paths(); GEDAI_alpha_figures('${STAGE}')"
