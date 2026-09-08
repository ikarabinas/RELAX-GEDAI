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

module load matlab

cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath('/home/nga4004/tools/eeglab2026.0.0/plugins/Fieldtrip-lite250523/'); \
addpath(pwd); addpath(fullfile(pwd,'tests')); \
run_integration_smoke"

LOG="logs/boundary_smoke_${SLURM_JOB_ID}.log"

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
