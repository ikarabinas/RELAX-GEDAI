% Save resting state EEG .mff files in .set format for use with RELAX preprocessing pipeline
% A standalone script for file conversion to EEGLAB format for when there are only a few files to process.
% Try run_set_conversion.py for batch processing of large datasets that need .mff -> .set conversion

N_CPU = 16;

% 
raw_dir = '/athena/grosenicklab/store/tms_eeg/mdd_dmpfc/subject32_m381_dmpfc_54/m381_dmpfc_day1';
save_dir= '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-dmpfc';

% Iterate over each resting state recording
if ~exist(save_dir, 'dir')
mkdir(save_dir);
end

% Initialize EEGLAB
eeglab;

raw_files = dir(fullfile(raw_dir, '**/*_reststate*.mff'));
for k = 1:numel(raw_files)
file = raw_files(k);
filepath = fullfile(file.folder, file.name);

EEG = pop_mffimport(filepath);
EEG = eeg_checkset(EEG);              % sanity check

% Convert data to double precision
%EEG.data = double(EEG.data);
%EEG = pop_editoptions(EEG, 'option_savetwofiles', 0, 'option_saveversion6', 1);

% Save as .set file with suffix
savefile = strrep(file.name, '.mff', '.set');
fprintf('%s', filepath)
%EEG = pop_saveset(EEG, 'filename', savefile, 'filepath', save_dir, 'savemode', 'twofiles');  % uncomment if you prefer .fdt + .set format
EEG = pop_saveset(EEG, 'filename', savefile, 'filepath', save_dir);  % save as one .set file by default

% Print the full save paths
full_savepath = fullfile(save_dir, savefile);
fprintf('Saving file to: %s\n', full_savepath);

end

eeglab redraw;

