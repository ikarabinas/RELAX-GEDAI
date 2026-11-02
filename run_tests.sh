#!/bin/bash
#SBATCH --job-name=relax-tests
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/tests_%j.log
# Unit tests, then the end-to-end smoke test. Submit from the repo or worktree root:
#   sbatch run_tests.sh
mkdir -p logs; module load matlab; cd "$SLURM_SUBMIT_DIR"
matlab -nodisplay -nosplash -batch "addpath('/home/nga4004/tools/eeglab2026.0.0'); addpath(pwd); addpath(fullfile(pwd,'tests')); eeglab nogui; failed = {}; try, run_all_tests; catch e, disp(getReport(e,'basic')); failed{end+1} = 'run_all_tests'; end; try, run_integration_smoke; disp('SMOKE PASS'); catch e, disp(getReport(e,'basic')); failed{end+1} = 'run_integration_smoke'; end; if ~isempty(failed), error('FAILED: %s', strjoin(failed, ', ')); end"
