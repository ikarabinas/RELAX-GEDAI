% GEDAI_alpha_masks() - Which samples the GEDAI alpha figures may count.
%
% Usage:
%   >> M = GEDAI_alpha_masks(EEG, opts);
%
% Inputs:
%   EEG  - EEGLAB continuous EEG structure; boundary and TMS_train_* events are read.
%   opts - struct from GEDAI_alpha_figures('opts'); uses edgeSec, minSegSec,
%          restTrainPadSec.
%
% Outputs (all masks 1 x pnts logical):
%   M.valid           - sample is more than edgeSec from any boundary or file edge, in a
%                       segment of at least minSegSec
%   M.rest            - valid, and before the first TMS_train_* event or after the last
%                       one plus restTrainPadSec
%   M.restIsWholeFile - true when the file has no TMS_train_* events
%   M.nExcluded       - number of samples that are not valid
%   M.nTinySegments   - segments shorter than 2 samples

function M = GEDAI_alpha_masks(EEG, opts)

    nPoints = EEG.pnts;
    srate = EEG.srate;
    segments = RELAX_get_segments(EEG);
    segmentLengths = segments(:, 2) - segments(:, 1) + 1;

    edgeSamples = round(opts.edgeSec * srate);
    minSegmentSamples = round(opts.minSegSec * srate);

    valid = false(1, nPoints);
    for s = find(segmentLengths >= minSegmentSamples)'
        valid(segments(s, 1) + edgeSamples : segments(s, 2) - edgeSamples) = true;
    end

    trainLatencies = [];
    for e = 1:numel(EEG.event)
        thisType = EEG.event(e).type;
        if (ischar(thisType) || isstring(thisType)) && startsWith(thisType, 'TMS_train_')
            trainLatencies(end+1) = EEG.event(e).latency; %#ok<AGROW>
        end
    end

    restIsWholeFile = isempty(trainLatencies);
    rest = valid;
    if ~restIsWholeFile
        sampleIndex = 1:nPoints;
        duringTrains = sampleIndex >= min(trainLatencies) & ...
            sampleIndex <= max(trainLatencies) + opts.restTrainPadSec * srate;
        rest = valid & ~duringTrains;
    end

    M.valid = valid;
    M.rest = rest;
    M.restIsWholeFile = restIsWholeFile;
    M.nExcluded = sum(~valid);
    M.nTinySegments = sum(segmentLengths < 2);
end
