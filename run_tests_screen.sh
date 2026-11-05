#!/bin/bash
#SBATCH --job-name=relax-tests-screen
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/tests_screen_%j.log
mkdir -p logs; module load matlab; cd "$SLURM_SUBMIT_DIR"
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); addpath(pwd); addpath(fullfile(pwd,'tests')); eeglab nogui; test_RELAX_screen_run_score; test_RELAX_screen_run_sheets; test_RELAX_screen_sheet_data; disp('SCREEN TESTS PASS')"
