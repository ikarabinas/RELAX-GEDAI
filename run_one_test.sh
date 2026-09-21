#!/bin/bash
#SBATCH --job-name=relax-one-test
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/one_test_%j.log
#SBATCH --error=logs/one_test_%j.error

# Run ONE tests/test_*.m function on a compute node.
#   mkdir -p logs && sbatch --wait run_one_test.sh test_RELAX_reference_channel
# `--wait` blocks until the job ends and returns its exit status, so a caller can chain
# `&& echo PASS`. The log is logs/one_test_<jobid>.log.

TEST_NAME="$1"
if [ -z "$TEST_NAME" ]; then echo "usage: sbatch --wait run_one_test.sh <test_function_name>" >&2; exit 2; fi

mkdir -p logs
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath(pwd); addpath(fullfile(pwd,'tests')); \
RELAX_add_vendored_noisetools(); \
${TEST_NAME}; disp('ONE_TEST_PASS')"
