% RELAX_pop_eye_mask() - Eye-movement periods on the raw timeline, for the pop stage.
%
% Usage:
%   >> [mask, info] = RELAX_pop_eye_mask(EEG, RELAX_cfg);
%
% The electrode-pop stage runs on the raw recording, long before the pipeline's own
% artifact marking. A horizontal saccade looks to the step detector like a local cluster
% of pops on the peri-ocular electrodes, and repairing it would remove EOG structure the
% later ICA relies on. This function builds, on a temporary copy, the same eye-movement
% masks the pipeline uses later - RELAX_blinks_IQR_method (blinks) and RELAX_horizontaleye
% (saccades) - so RELAX_flag_eye_steps can leave those clusters alone.
%
% Both detectors are given what they would see in the pipeline: the copy is average
% re-referenced and band-passed 1-25 Hz (the blink method filters its own electrodes; the
% HEOG electrodes are filtered here). Only the mask is kept: the input EEG is not
% modified, no events are added to it, and every failure (missing electrodes, detector
% error) degrades to an all-false contribution with a warning, so at worst a saccade
% cluster gets repaired like a pop.
%
% The HEOG contribution is OFF unless RELAX_cfg.PopEye_IncludeHEOG is set (since
% 2026-09-18). RELAX_horizontaleye looks at two channels only and asks for opposite-sign
% excursions; after the zero-phase 1-25 Hz bandpass a DC step is a bipolar excursion
% (negative lobe before, positive after), so two same-sign pops a few tens of ms apart on
% the HEOG pair satisfy that test and mark +-HorizontalEyeMovementFocus around themselves -
% the pops generate the very mask that then exempts them from repair. The pair is still
% resolved into info.heogChannels either way, because the wrapper counts it among the eye
% electrodes that are protected from pop-based rejection.
%
% Inputs:
%   EEG       - raw continuous EEGLAB struct (native sample rate, before any filter)
%   RELAX_cfg - needs BlinkElectrodes, HEOGLeftpattern / HEOGRightpattern and the usual
%               blink / horizontal-eye fields; missing ones take RELAX's defaults.
%               PopEye_IncludeHEOG (0/1, default 0) adds the horizontal-eye detector.
%
% Outputs:
%   mask - logical [1 x EEG.pnts], true during blink (and, if enabled, horizontal-eye-
%          movement) periods
%   info - .blinkStatus / .heogStatus ('ok', 'disabled', or the reason that detector
%          contributed nothing), .blinkFraction / .heogFraction (fraction of samples each
%          marked), .maskFraction, .heogChannels (the resolved HEOG pair, {} if not
%          present in the montage), .rereferenced

function [mask, info] = RELAX_pop_eye_mask(EEG, RELAX_cfg)

    pnts = EEG.pnts;
    blinkMask = false(1, pnts);
    heogMask  = false(1, pnts);
    info = struct('blinkStatus', '', 'heogStatus', '', 'blinkFraction', 0, ...
                  'heogFraction', 0, 'maskFraction', 0, 'heogChannels', {{}}, ...
                  'rereferenced', false);

    if ~isfield(RELAX_cfg, 'ms_per_sample'); RELAX_cfg.ms_per_sample = 1000 / EEG.srate; end
    if ~isfield(RELAX_cfg, 'HighPassFilter'); RELAX_cfg.HighPassFilter = 1;   end
    if ~isfield(RELAX_cfg, 'LowPassFilter');  RELAX_cfg.LowPassFilter  = 200; end
    if ~isfield(RELAX_cfg, 'RespectBoundaries'); RELAX_cfg.RespectBoundaries = 0; end
    if ~isfield(RELAX_cfg, 'PopEye_IncludeHEOG'); RELAX_cfg.PopEye_IncludeHEOG = 0; end
    labels = {EEG.chanlocs.labels};

    % Temporary copy, average re-referenced the way RELAX_blinks_IQR_method would do it
    % itself. ListOfChannelsAfterRejections / allchan are set to the channels present so
    % RELAX_average_rereference neither interpolates anything nor drops anything.
    tmp = EEG;
    tmp.RELAX = struct('ListOfChannelsAfterRejections', {labels});
    tmp.allchan = EEG.chanlocs;
    tmp.RELAXProcessing = struct();
    try
        tmp = evalc_quiet(@() RELAX_average_rereference(tmp));
        info.rereferenced = true;
    catch err
        warning('RELAX_pop_eye_mask:reref', ...
            'average re-reference of the temporary copy failed (%s); detectors run on the raw reference', err.message);
        tmp = EEG;
        tmp.RELAX = struct('ListOfChannelsAfterRejections', {labels}, 'Data_has_been_averagerereferenced', 1);
        tmp.allchan = EEG.chanlocs;
        tmp.RELAXProcessing = struct();
    end

    % --- blinks --------------------------------------------------------------------------
    if ~isfield(RELAX_cfg, 'BlinkElectrodes') || sum(ismember(RELAX_cfg.BlinkElectrodes(:), labels)) < 2
        info.blinkStatus = 'fewer than two blink electrodes present';
    else
        try
            b = evalc_quiet(@() RELAX_blinks_IQR_method(tmp, tmp, RELAX_cfg));
            if isfield(b, 'RELAX') && isfield(b.RELAX, 'eyeblinkmask')
                blinkMask = reshape(b.RELAX.eyeblinkmask == 1, 1, []);
                info.blinkStatus = 'ok';
            else
                info.blinkStatus = 'blink method returned no mask';
            end
        catch err
            info.blinkStatus = err.message;
        end
    end
    if ~strcmp(info.blinkStatus, 'ok')
        warning('RELAX_pop_eye_mask:blinks', 'no blink mask for the pop stage: %s', info.blinkStatus);
    end

    % --- horizontal eye movements -------------------------------------------------------
    left  = firstPresent(RELAX_cfg, 'HEOGLeftpattern',  labels);
    right = firstPresent(RELAX_cfg, 'HEOGRightpattern', labels);
    if isempty(left) || isempty(right)
        info.heogStatus = 'no lateral HEOG electrode pair present';
    else
        info.heogChannels = {left, right};
        if ~RELAX_cfg.PopEye_IncludeHEOG
            info.heogStatus = 'disabled';
        else
            try
                h = evalc_quiet(@() pop_select(tmp, 'channel', {left, right}));
                h = evalc_quiet(@() RELAX_filtbutter(h, 1, 25, 4, 'bandpass', 'acausal', RELAX_cfg.RespectBoundaries));
                h.RELAXProcessing.Details.NaNsForNonEvents    = zeros(1, pnts);
                h.RELAXProcessing.Details.NoiseMaskFullLength = zeros(1, pnts);
                h = evalc_quiet(@() RELAX_horizontaleye(h, RELAX_cfg));
                heogMask = reshape(h.RELAXProcessing.Details.horizontaleyemask == 1, 1, []);
                info.heogStatus = 'ok';
            catch err
                info.heogStatus = err.message;
            end
        end
    end
    if ~any(strcmp(info.heogStatus, {'ok', 'disabled'}))
        warning('RELAX_pop_eye_mask:heog', 'no horizontal-eye mask for the pop stage: %s', info.heogStatus);
    end

    mask = blinkMask | heogMask;
    info.blinkFraction = mean(blinkMask);
    info.heogFraction  = mean(heogMask);
    info.maskFraction  = mean(mask);
end

function lab = firstPresent(cfg, field, labels)
% First label of a priority pattern that exists in the montage ('' if none), matching
% RELAX_horizontaleye's own case-insensitive search order.
    lab = '';
    if ~isfield(cfg, field); return; end
    pattern = cellstr(cfg.(field));
    for k = 1:numel(pattern)
        hit = find(strcmpi(labels, pattern{k}), 1);
        if ~isempty(hit); lab = labels{hit}; return; end
    end
end

function out = evalc_quiet(fn) %#ok<INUSD> (fn is used inside the evalc string)
% Run fn() with its command-window chatter captured (the detectors print several lines
% per call that mean nothing at this stage).
    [~, out] = evalc('fn()');
end
