#!/bin/bash
#SBATCH --job-name=relax-run-all-tests
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/run_all_tests_%j.log
#SBATCH --error=logs/run_all_tests_%j.error

# Full unit-test suite (tests/test_*.m via tests/run_all_tests.m).
# Submit with: mkdir -p logs && sbatch run_all_tests.sh

mkdir -p logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath(pwd); addpath(fullfile(pwd,'tests')); \
run_all_tests"
