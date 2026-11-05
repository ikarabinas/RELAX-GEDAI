#!/bin/bash
#SBATCH --job-name=extreme-survey
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=1-89%30
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.error

# One Extremes_Rejected file per array task (index into the sorted file list).
# Submit from this directory:  mkdir -p logs && sbatch run_survey.sh
# Then:                        sbatch --dependency=afterany:<jobid> run_summarise.sh
# K_OVERRIDE=20 RESULTS=results_k20 reruns RELAX with that k instead of the run's.

set -eo pipefail
RUN_DIR="${RUN_DIR:-/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26_2exp/rest_treat_joined/RELAXProcessed_gentlerICA_picard}"
here=$(pwd)
repo=$(cd .. && pwd)
module load matlab

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0/'); eeglab nogui; \
addpath('/home/nga4004/tools/eeglab2026.0.0/plugins/Fieldtrip-lite250523/'); \
addpath('$repo'); addpath('$here'); \
c = parcluster('local'); jsl = fullfile(getenv('TMPDIR'), 'mjs_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}'); \
mkdir(jsl); c.JobStorageLocation = jsl; parpool(c, ${SLURM_CPUS_PER_TASK}); \
survey_extreme_rejection_file('$RUN_DIR', ${SLURM_ARRAY_TASK_ID}, '$here/${RESULTS:-results}', [${K_OVERRIDE:-}])"
