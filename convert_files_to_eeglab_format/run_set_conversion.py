import os
import pandas as pd
import subprocess

'''
Convert all .mff EEG files for a given participant list to .set files for use with EEGLAB
Here, the data have been stored in a directory structure as follows: .../tms_eeg/<tms_target>/<ppt_id>/<treatment_day>/<mff_file>
Remember to load MATLAB in the terminal with "module load matlab/R2023a" or an equivalent command before running this script with "python run_set_conversion.py"
'''
# Define data paths and read in participant IDs csv
DATA_DIR="/athena/grosenicklab/store/tms_eeg/ocd_lofc/"
SAVE_DIR="/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-ofc/"
ppts_csv_path="/home/imk2003/Documents/updated_subject_list_lofc.csv"
matlab_script = '/home/imk2003/Documents/MATLAB/eeglab/plugins/RELAX/convert_files_to_eeglab_format/save_mff_to_set.m'
ppts_csv = pd.read_csv(ppts_csv_path)
days_list = ['day1', 'day2', 'day3', 'day4', 'day5', 'week', 'baseline']

# Extract TMS target info from data dir
data_dir_name = os.path.basename(os.path.normpath(DATA_DIR))  # "mdd_dlpfc"
tms_target = data_dir_name.split('_')[-1].upper()             # "DLPFC"

# Filter CSV to keep only rows with specified TMS target and extract matching ppts
#filtered_csv = ppts_csv[ppts_csv['tms_target'].str.upper() == tms_target]
filtered_csv = ppts_csv
ppts_list = list(filtered_csv['record_id'])

# Extract unique participant IDs from saved files
saved_files = os.listdir(SAVE_DIR)
already_converted_ppt_files = set()
for file in saved_files:
    if file.endswith('.set'):
        if 'crossover' in file:
            ppt_id = file.split('_')[0] + '_crossover'
        else:
            ppt_id = file.split('_')[0]
            
        already_converted_ppt_files.add(ppt_id)


# Run save_mff_to_set for each ppt-day resting state EEG file path
def run_saveset(ppt_day_datapath):
    # MATLAB command string (call with two arguments)
    matlab_cmd = f"save_mff_to_set('{ppt_day_datapath}', '{SAVE_DIR}')"
    try:
        subprocess.run(['matlab', '-batch', matlab_cmd], check=True)
    except subprocess.CalledProcessError as e:
        print(f'[ERROR] Error during preprocessing: {e}')      


# Takes a participant record id and returns the file path corresponding to the specified analysis day dir
def get_eeg_day_path(ppt_id, analysis_day):

    # Filter data directories for specified participant and analysis day
    try:
        ppt_dir = str([dir for dir in os.listdir(DATA_DIR) if ppt_id in dir][0])
        ppt_dirpath = os.path.join(DATA_DIR, ppt_dir)

        day_path = str([dir for dir in os.listdir(ppt_dirpath) if analysis_day in dir][0])
        ppt_day_dirpath = os.path.join(ppt_dirpath, day_path)
    except IndexError:
        if 'crossover' in ppt_id:
            stripped_id = ppt_id.split('_')[0]  # e.g. 'c110_crossover' -> 'c110'
            try:
                ppt_dir = str([dir for dir in os.listdir(DATA_DIR) if stripped_id in dir][0])
                ppt_dirpath = os.path.join(DATA_DIR, ppt_dir)
                day_path = str([dir for dir in os.listdir(ppt_dirpath) if analysis_day in dir][0])
                ppt_day_dirpath = os.path.join(ppt_dirpath, day_path)
            except IndexError:
                raise FileNotFoundError(f'No data found for {ppt_id} {analysis_day} (also tried {stripped_id}).')
        else:
            raise FileNotFoundError(f'No data found for {ppt_id} {analysis_day}.')
    return ppt_day_dirpath

# Run .set file conversion for each participant
for ppt in ppts_list:
    if ppt in already_converted_ppt_files:  # skip already processed files
        print(f"Skipping {ppt}: already processed.")
        continue
    else:
        for day in days_list:
            try:
                ppt_day_datapath = get_eeg_day_path(ppt, day)
            # Catch when files do not exist
            except FileNotFoundError as e:
                print(f"[WARNING] {e}. Skipping {ppt} {day}.")
                continue   # move on to next ppt-day file

            print(ppt_day_datapath)
            run_saveset(ppt_day_datapath)
