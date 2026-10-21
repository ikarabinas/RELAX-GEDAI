function [L, ch_names] = RELAX_load_ppt_leadfield(leadfield_dir, eeg_filepath)
% RELAX_LOAD_PPT_LEADFIELD Load the participant leadfield matching a recording filename.
%
% Finds the leadfield .mat in leadfield_dir named '<participant>_*.mat', where
% <participant> is the first underscore-delimited part of the recording filename
% (e.g. 'm155' for 'm155_dlpfc_day3_stitched.set'). These files are produced by
% generate_leadfield.py in ~/TMS_response_variability and hold the raw
% (non-average-referenced) leadfield in 'L' (nchan x ndipoles) and the channel labels
% in 'ch_names'. Copied from load_ppt_leadfield in ~/TMS_response_variability/iTBS_pipeline.m.
%
% Output:
%   L        - raw leadfield, nchan x ndipoles (double)
%   ch_names - 1 x nchan cell array of channel labels corresponding to rows of L
    [~, fname_noext, ~] = fileparts(eeg_filepath);
    parts = strsplit(fname_noext, '_');
    prefix = [parts{1} '_'];

    dirList = dir(fullfile(leadfield_dir, '*.mat'));
    if isempty(dirList)
        error('RELAX_load_ppt_leadfield:noFiles', 'No .mat files found in %s', leadfield_dir);
    end
    matchIdx = find(startsWith({dirList.name}, prefix));
    if isempty(matchIdx)
        error('RELAX_load_ppt_leadfield:noMatch', 'No matching leadfield found for prefix: %s', prefix);
    elseif numel(matchIdx) > 1
        error('RELAX_load_ppt_leadfield:multipleMatches', 'Multiple matching leadfield files found for prefix: %s', prefix);
    end

    lf_path = fullfile(leadfield_dir, dirList(matchIdx).name);
    fprintf('Loading participant leadfield from: %s\n', lf_path);
    data = load(lf_path);
    L = double(data.L);
    ch_names = data.ch_names;
    if ~iscell(ch_names)
        ch_names = cellstr(ch_names);
    end
    ch_names = ch_names(:)';
end
