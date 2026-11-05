#!/bin/bash
#SBATCH --job-name=extreme-summary
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.error

set -eo pipefail
here=$(pwd)
module load matlab
matlab -nodisplay -nosplash -batch "addpath('$here'); summarise_extreme_rejection('$here/${RESULTS:-results}')" \
    | tee "$here/summary_${RESULTS:-results}.txt"
