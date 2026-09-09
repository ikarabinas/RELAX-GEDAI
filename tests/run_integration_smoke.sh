#!/bin/bash
#SBATCH --job-name=relax-boundary-smoke
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/boundary_smoke_%j.log
#SBATCH --error=logs/boundary_smoke_%j.error

# End-to-end smoke run of the boundary-aware pipeline on the cropped fixture.
# Submit with: sbatch tests/run_integration_smoke.sh
#
# NOTE: the #SBATCH --output/--error directives above are evaluated by sbatch at
# submission time, BEFORE any line of this script runs - so logs/ must already exist
# at submission time or sbatch fails outright. The `mkdir -p logs` below only helps a
# script that is re-run after logs/ has been removed; it cannot help the directives
# above on first submission from a fresh clone. Run `mkdir -p logs` once yourself
# (or `cd` to a checkout where it already exists) before `sbatch tests/run_integration_smoke.sh`.

mkdir -p logs

module load matlab

cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath('/home/nga4004/tools/eeglab2026.0.0/plugins/Fieldtrip-lite250523/'); \
addpath(pwd); addpath(fullfile(pwd,'tests')); \
run_integration_smoke"
MATLAB_STATUS=$?

LOG="logs/boundary_smoke_${SLURM_JOB_ID}.log"

# Fail on MATLAB's own exit status FIRST. The greps below only prove PREP was skipped,
# which is logged early in the run - they say nothing about whether the rest of the
# pipeline or the assert_smoke_outcome checks succeeded.
if [ $MATLAB_STATUS -ne 0 ]; then
    echo "FAIL: MATLAB exited with status $MATLAB_STATUS" >&2
    exit $MATLAB_STATUS
fi

# The in-MATLAB assertions cannot distinguish "PREP was skipped" from "PREP ran and found
# nothing" - both leave PREPBasedChannelToReject empty. The log can.
if ! grep -q "Skipping PREP findNoisyChannels" "$LOG"; then
    echo "FAIL: log does not show PREP being skipped" >&2
    exit 1
fi
if grep -q "Removing noisy channels:" "$LOG"; then
    echo "FAIL: log shows PREP rejected channels - it was not skipped" >&2
    exit 1
fi
echo "Confirmed from log: PREP findNoisyChannels was skipped."
