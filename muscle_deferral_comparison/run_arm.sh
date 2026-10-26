#!/bin/bash
#SBATCH --job-name=muscle-defer
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.error

# One arm of the muscle-screen placement comparison, one participant per array task
# (indices from <myPath>/participant_list.mat). Submit from this directory:
#   sbatch --array=2,5,6,10,15,16,20 --export=ALL,ARM=Before run_arm.sh
#   sbatch --array=2,5,6,10,15,16,20 --export=ALL,ARM=After  run_arm.sh
# Before = RejectMuscleChannelsAfterCleaning 0 (RELAX default), After = 1. Everything
# else is the branch's RELAX_SET_PARAMETERS_AND_RUN_slurm.m, run against this worktree's
# code rather than the main checkout's, writing to <myPath>/RELAXProcessed_muscle<ARM>.

set -eo pipefail
here=$(pwd)
worktree=$(cd .. && pwd)
case "$ARM" in
    Before) deferred=0 ;;
    After)  deferred=1 ;;
    *) echo "ARM must be Before or After, got '$ARM'" >&2; exit 1 ;;
esac

mkdir -p arms "logs/$ARM"
entry="arms/RELAX_SET_PARAMETERS_AND_RUN_muscle${ARM}_${SLURM_ARRAY_TASK_ID}.m"
sed -e "s#addpath('/home/nga4004/gitRepos/RELAX-GEDAI/');#addpath('$worktree/');#" \
    -e "s#RELAX_cfg.OutputFolderName='RELAXProcessed_broadbandGEDAI';#RELAX_cfg.OutputFolderName='RELAXProcessed_muscle$ARM';#" \
    -e "s#RELAX_cfg.RejectMuscleChannelsAfterCleaning=0;#RELAX_cfg.RejectMuscleChannelsAfterCleaning=$deferred;#" \
    -e "s#sprintf('%s/slurm_task_%d.log', RELAX_cfg.myPath, task_id)#sprintf('%s/slurm_task_%d.log', '$here/logs/$ARM', task_id)#" \
    ../RELAX_SET_PARAMETERS_AND_RUN_slurm.m > "$entry"

# Every substitution must have landed, or the arm would silently run the main checkout's
# code, write over the broadbandGEDAI outputs, or run the wrong placement.
grep -qF "addpath('$worktree/');" "$entry"
grep -qF "RELAX_cfg.OutputFolderName='RELAXProcessed_muscle$ARM';" "$entry"
grep -qF "RELAX_cfg.RejectMuscleChannelsAfterCleaning=$deferred;" "$entry"
[ "$(grep -cF "$here/logs/$ARM" "$entry")" -eq 2 ]
if grep -qF "gitRepos/RELAX-GEDAI/');" "$entry"; then exit 1; fi

module load matlab/R2023a
matlab -nodisplay -nosplash -nodesktop -batch "run('$here/$entry')"
