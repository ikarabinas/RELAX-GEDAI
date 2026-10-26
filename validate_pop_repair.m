function results = validate_pop_repair(outFile, opts)
% VALIDATE_POP_REPAIR - Measure what the two-window electrode-pop stage does on real files.
%
% Reproduces RELAX_Wrapper_slurm's pop stage (detector -> bad-contact rejection -> eye
% mask -> spatial fill + common-mode excision) on the four surveyed recordings, on the
% production 204-channel montage, and reports for each file:
%
%   1. detection load: steps, events, common-mode / eye-movement / rejected-channel
%      counts, interpolated fraction per channel (ceiling 2 %), excised seconds;
%   2. what the pipeline's 1 Hz highpass leaves at each filled transition (std of the
%      highpassed channel over [onset - 1 s, offset + 1 s]) for four versions of the
%      data on a shared timeline (excision off): RAW (pop left in), LEGACY (the previous
%      derivative detector + single-sample-cut repair at commit f6649a8, evaluated at the
%      new detector's sites), LINEAR bridge, SPATIAL fill; each against a step-free
%      control window on the same channel;
%   3. a leave-one-out check of the spatial fill where the truth is known: on a step-free
%      control window of the same channel and length as each fill, is the neighbours'
%      spline estimate closer to the channel's actual anchor-relative detail than the
%      linear chord (zero detail) is?
%   4. on one file, ICA agreement raw vs repaired and raw vs detrend-only (the ladder with
%      no fills, to separate the two), plus split-half reproducibility (first vs second
%      half of the recording, decomposed separately) for raw and repaired - a reference
%      that needs no pop-free data, which a poppy file does not have (masking +/- 1 s
%      around every fill removed 73 % of m159): 1 Hz highpass, 250 Hz, PICARD on the same
%      channels,
%      ICs matched by |corr| of the mixing columns (matchpairs), overall and per ICLabel
%      class of the raw IC, ICLabel class histogram per condition; plus the floor of the
%      statistic, raw against raw minus its last 1 % of samples.
%
% Run on a compute node, never a login node:
%   sbatch run_validate_pop_repair.sh
% or interactively inside an allocation (srun does not start a login shell, so the module
% load has to happen inside it):
%   srun -p scu-cpu -c 4 --mem=48G --time=03:00:00 bash -c '
%     module load matlab; cd ~/gitRepos/RELAX-GEDAI
%     matlab -nodisplay -nosplash -batch "addpath(\"~/tools/eeglab2026.0.0\"); addpath(pwd); validate_pop_repair(\"pop_validation.mat\")"'
%
% opts (all optional): .files, .dataDir, .icaFile (index into files; 0 = skip ICA),
%   .icaPCA (80), .legacy (true), .montageFile, .refineTransition (true, matching the
%   shipped PopDetect_RefineTransition default; pass false for the unrefined-transition arm).

    if nargin < 1 || isempty(outFile); outFile = 'pop_validation.mat'; end
    if nargin < 2; opts = struct(); end
    opts = setDefault(opts, 'dataDir', '~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Sep26/rest_treat_joined');
    opts = setDefault(opts, 'files', { ...
        'm159_dlpfc_day1_stitched.set', ...
        'm497_dlpfc_day1_stitched.set', ...
        'm208_dlpfc_day1_stitched.set', ...
        'm583_dlpfc_day1_stitched.set'});
    opts = setDefault(opts, 'montageFile', '/home/nga4004/TMS_response_variability/montage_204_egi_indices.txt');
    opts = setDefault(opts, 'icaFile', 1);
    opts = setDefault(opts, 'icaPCA', 80);
    opts = setDefault(opts, 'legacy', true);
    opts = setDefault(opts, 'refineTransition', true);

    eeglab nogui;
    RELAX_add_vendored_noisetools();
    cfg = stageConfig(opts.montageFile);
    cfg.PopDetect_RefineTransition = opts.refineTransition;
    legacy = legacyTools(opts.legacy);

    results = {};      % one struct per file; fields differ (ICA / legacy arms), so a cell
    for k = 1:numel(opts.files)
        fprintf('\n===== %s\n', opts.files{k});
        EEG = pop_loadset('filename', opts.files{k}, 'filepath', opts.dataDir);
        EEG = pop_select(EEG, 'nochannel', intersect(cfg.ElectrodesToDelete, {EEG.chanlocs.labels}));
        EEG = eeg_checkset(EEG);
        EEG.allchan = EEG.chanlocs;
        minutes = EEG.pnts / EEG.srate / 60;

        r = struct('file', opts.files{k}, 'nbchan', EEG.nbchan, 'minutes', minutes);

        % ---- 1. the production path -------------------------------------------------
        tStage = tic;
        S = runStage(EEG, cfg);
        r.stageSeconds = toc(tStage);
        r.detect = S.detect; r.repair = S.repair; r.eye = S.eyeInfo; r.rejected = S.rejected;
        r.steps  = S.steps;
        printStage(r, S);

        % ---- 2. residual after the 1 Hz highpass, four versions, shared timeline ------
        base = S.preRepairEEG;                                   % rejected channels gone
        fillList = fillIntervals(S.steps);                       % [ch lo hi] on base
        V = struct();
        V.raw     = base;
        V.linear  = repairWith(base, S.steps, S.basisExclude, cfg, 'linear');
        V.spatial = repairWith(base, S.steps, S.basisExclude, cfg, 'spatial');
        if ~isempty(legacy)
            try
                lsteps = legacy.detect(base);
                V.legacy = legacy.repair(base, lsteps);
                r.legacyNSteps = numel(lsteps.sample);
            catch err
                warning('legacy stage failed on %s: %s', opts.files{k}, err.message);
            end
        end
        r.ringing = ringingAtSites(V, fillList, base.srate);
        printRinging(r.ringing);

        % ---- 3. leave-one-out gain of the spatial estimate on step-free windows -------
        r.loo = spatialLeaveOneOut(V.spatial, fillList, S.basisExclude, cfg.PopRepair_AnchorMs);
        fprintf(['LOO on %d step-free control windows: rms detail %.2f uV (linear residual) -> ' ...
                 '%.2f uV after spline estimate; spatial <= linear on %.0f%% (%.0f%% of windows > 50 ms)\n'], ...
            r.loo.n, r.loo.rmsLinear, r.loo.rmsSpatial, 100*r.loo.fracSpatialBetter, 100*r.loo.fracSpatialBetterLong);

        % ---- 4. ICA agreement on one file ---------------------------------------------
        if k == opts.icaFile
            fprintf('ICA comparison (PICARD, pca %d) raw vs repaired vs detrend-only ...\n', opts.icaPCA);
            try
                noSteps = structfun(@(f) f(false(size(f, 1), 1), :), S.steps, 'UniformOutput', false);
                I = struct('raw', base, 'repaired', S.EEG, ...
                           'detrendOnly', repairWith(base, noSteps, S.basisExclude, cfg, 'linear'));
                r.ica = compareICA(I, opts.icaPCA);
                for nm = {'repaired', 'detrendOnly', 'floor'}
                    c = r.ica.(nm{1});
                    fprintf('ICA raw vs %-11s: %d matched pairs, median |r| %.3f, %d with |r| > 0.9, %d with |r| > 0.95\n', ...
                        nm{1}, c.nPairs, c.medianAbsR, c.nAbove090, c.nAbove095);
                end
                fprintf('ICA raw vs repaired by raw ICLabel class: %s\nICLabel classes raw / repaired: %s\n', ...
                    r.ica.byClassSummary, r.ica.classSummary);
                fprintf('ICA pairwise median |r| (rows/cols %s):\n', strjoin(r.ica.arms, ' '));
                disp(round(r.ica.pairMedianAbsR, 3));
                for nm = {'raw', 'repaired'}
                    c = r.ica.splitHalf.(nm{1});
                    fprintf('ICA split-half (first vs second half) %-8s: median |r| %.3f, %d with |r| > 0.9, %d with |r| > 0.95\n', ...
                        nm{1}, c.medianAbsR, c.nAbove090, c.nAbove095);
                end
            catch err
                warning('validate_pop_repair:ica', 'ICA comparison failed: %s', err.message);
                r.ica = struct('error', err.message);
            end
        end

        results{end+1} = r; %#ok<AGROW>
        save(outFile, 'results', 'opts', 'cfg');
    end

    printSummary(results);
    save(outFile, 'results', 'opts', 'cfg');
    fprintf('\nSaved %s\n', outFile);
end

% ==========================================================================================
function cfg = stageConfig(montageFile)
% The subset of RELAX_SET_PARAMETERS_AND_RUN_slurm.m the pop stage reads.
    cfg = struct();
    cfg.RespectBoundaries = 1;
    cfg.Perform_pop_artifact_repair = 1;
    cfg.BlinkElectrodes = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
    cfg.HEOGLeftpattern  = ["E46", "E54", "E252", "E47", "E248", "E241", "E244"];
    cfg.HEOGRightpattern = ["E10", "E1", "E226", "E2", "E230", "E238", "E234"];
    cfg.HorizontalEyeMovementType = 2;
    cfg.HorizontalEyeMovementThreshold = 2;
    cfg.HorizontalEyeMovementThresholdIQR = 1.5;
    cfg.HorizontalEyeMovementTimepointsExceedingThreshold = 25;
    cfg.HorizontalEyeMovementTimepointsTestWindow = 49;
    cfg.HorizontalEyeMovementFocus = 200;
    cfg.HighPassFilter = 1; cfg.LowPassFilter = 200;
    cfg.Perform_impedance_bad_channel_rejection = 1;
    cfg.MaxImpedanceKOhm = 500;
    cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.10;
    cfg.ElectrodesToDelete = {'E257'};
    if ~isempty(montageFile) && exist(montageFile, 'file') == 2
        idx = readmatrix(montageFile); idx = idx(~isnan(idx));
        drop = setdiff(1:257, idx(:)');
        cfg.ElectrodesToDelete = [arrayfun(@(n) sprintf('E%d', n), drop(:), 'UniformOutput', false); cfg.ElectrodesToDelete(:)];
    end
    % stage parameters: the wrapper's defaults (kept in one place there)
    cfg.PopDetect_WindowMs = 101; cfg.PopDetect_TransitionMs = 6;
    cfg.PopDetect_ZThreshold = 10; cfg.PopDetect_ChannelZThreshold = 20;
    cfg.PopDetect_MinAmplitudeMicrovolts = 15; cfg.PopDetect_RefractoryMs = 10;
    cfg.PopDetect_RefineTransition = true;
    cfg.PopDetect_StepsPerMinuteToRejectChannel = 6;
    cfg.PopDetect_MaxTransitionFractionToRejectChannel = 0.02;
    cfg.PopDetect_CommonModeChannelFraction = 0.15;
    cfg.PopDetect_SlopeBaselineMs = 500; cfg.PopDetect_SlopeNoiseFloorMicrovolts = 1;
    cfg.PopDetect_SlopeSpreadQuantile = 0.95;
    cfg.PopEye_MaxAmplitudeMicrovolts = 500; cfg.PopEye_IncludeHEOG = 0;
    cfg.PopRepair_MaxPolyOrder = 3; cfg.PopRepair_MaxCutoffHz = 1.0; cfg.PopRepair_MaxWindowSeconds = 2.0;
    cfg.PopRepair_DetrendThresh = 3; cfg.PopRepair_DetrendIterations = 3;
    cfg.PopRepair_AnchorMs = 3; cfg.PopRepair_FillMethod = 'spatial'; cfg.PopRepair_MinBasisChannels = 16;
    cfg.PopRepair_ExciseCommonMode = false; cfg.PopRepair_MaxExciseMs = 500; cfg.PopRepair_MinSegmentSeconds = 1.0;
end

function S = runStage(EEG, cfg)
% The wrapper's pop stage (RELAX_Wrapper_slurm.m, "Detect and repair electrode pops"),
% returning the intermediates the measurements need. Two wrapper steps are not reproduced.
% The impedance screen does not run, so high-impedance channels stay in the montage the
% detector sees (they are only kept out of the spline basis) and the channel budget starts
% unspent. The reference-pop pass does not run either.
    cfg.ms_per_sample = 1000 / EEG.srate;
    EEG = RELAX_channel_budget(EEG, 'init', cfg.MaxProportionOfElectrodesThatCanBeDeleted);
    [steps, dInfo] = RELAX_detect_steps(EEG, ...
        'WindowMs', cfg.PopDetect_WindowMs, 'TransitionMs', cfg.PopDetect_TransitionMs, ...
        'ZThreshold', cfg.PopDetect_ZThreshold, 'ChannelZThreshold', cfg.PopDetect_ChannelZThreshold, ...
        'MinAmplitudeMicrovolts', cfg.PopDetect_MinAmplitudeMicrovolts, 'RefractoryMs', cfg.PopDetect_RefractoryMs, ...
        'RefineTransition', cfg.PopDetect_RefineTransition, ...
        'StepsPerMinuteToRejectChannel', cfg.PopDetect_StepsPerMinuteToRejectChannel, ...
        'MaxTransitionFractionToRejectChannel', cfg.PopDetect_MaxTransitionFractionToRejectChannel, ...
        'CommonModeChannelFraction', cfg.PopDetect_CommonModeChannelFraction, ...
        'SlopeBaselineMs', cfg.PopDetect_SlopeBaselineMs, ...
        'SlopeNoiseFloorMicrovolts', cfg.PopDetect_SlopeNoiseFloorMicrovolts, ...
        'SlopeSpreadQuantile', cfg.PopDetect_SlopeSpreadQuantile);
    if ~dInfo.available; error('detector unavailable: %s', dInfo.status); end

    % eye movements first (on the complete montage), then bad contacts ranked over the
    % remaining steps - same order as the wrapper
    [mask, maskInfo] = RELAX_pop_eye_mask(EEG, cfg);
    eyeLabels = union(reshape(cfg.BlinkElectrodes, [], 1), reshape(maskInfo.heogChannels, [], 1));
    [steps, eyeInfo] = RELAX_flag_eye_steps(steps, mask, eyeLabels, ...
        'MaxAmplitudeMicrovolts', cfg.PopEye_MaxAmplitudeMicrovolts);
    eyeInfo.mask = maskInfo; eyeInfo.labels = eyeLabels;
    dInfo.badChannelsBeforeEyeFlag = dInfo.badChannels;
    dInfo = RELAX_rank_pop_channels(steps, dInfo);
    [steps, dInfo] = RELAX_protect_eye_channels(steps, dInfo, eyeLabels);

    [rejected, EEG] = RELAX_channel_budget(EEG, 'spend', 'pop', dInfo.badChannels);
    if ~isempty(rejected)
        keep = ~ismember(steps.label, rejected);
        steps = structfun(@(f) f(keep), steps, 'UniformOutput', false);
        EEG = pop_select(EEG, 'nochannel', rejected); EEG = eeg_checkset(EEG);
        [~, steps.channel] = ismember(steps.label, {EEG.chanlocs.labels});
    end

    basisExclude = reshape(dInfo.degenerateChannels, 1, []);
    if isfield(EEG.chanlocs, 'impedance')
        imp = NaN(1, EEG.nbchan);
        for c = 1:EEG.nbchan
            v = EEG.chanlocs(c).impedance;
            if isnumeric(v) && isscalar(v); imp(c) = double(v); end
        end
        basisExclude = union(basisExclude, {EEG.chanlocs(imp > cfg.MaxImpedanceKOhm).labels});
    end

    S = struct('detect', dInfo, 'steps', steps, 'rejected', {rejected}, 'eyeInfo', eyeInfo, ...
               'basisExclude', {basisExclude}, 'preRepairEEG', EEG);
    [S.EEG, S.repair] = repairWith(EEG, steps, basisExclude, cfg, cfg.PopRepair_FillMethod, cfg.PopRepair_ExciseCommonMode);
    if ~S.repair.available; error('repair unavailable: %s', S.repair.status); end
end

function [EEG, info] = repairWith(EEG, steps, basisExclude, cfg, fillMethod, excise)
    if nargin < 6; excise = false; end
    [EEG, info] = RELAX_repair_steps(EEG, steps, ...
        'MaxPolyOrder', cfg.PopRepair_MaxPolyOrder, 'MaxCutoffHz', cfg.PopRepair_MaxCutoffHz, ...
        'MaxWindowSeconds', cfg.PopRepair_MaxWindowSeconds, 'DetrendThresh', cfg.PopRepair_DetrendThresh, ...
        'DetrendIterations', cfg.PopRepair_DetrendIterations, 'AnchorMs', cfg.PopRepair_AnchorMs, ...
        'FillMethod', fillMethod, 'BasisExcludeLabels', basisExclude, ...
        'MinBasisChannels', cfg.PopRepair_MinBasisChannels, ...
        'ExciseCommonMode', excise, 'MaxExciseMs', cfg.PopRepair_MaxExciseMs, ...
        'MinSegmentSeconds', cfg.PopRepair_MinSegmentSeconds);
end

function F = fillIntervals(steps)
% [ch lo hi] of the steps the repair fills (not common-mode, not eye-movement).
    keep = ~steps.commonMode;
    if isfield(steps, 'eyeMovement'); keep = keep & ~steps.eyeMovement; end
    F = [steps.channel(keep), steps.onset(keep), steps.offset(keep)];
end

function printStage(r, S)
    d = r.detect; p = r.repair; e = r.eye;
    fprintf('stage %.0f s: %d steps on %d/%d channels (%.1f/min), %d events; median rise %.1f ms; median |A| %.0f uV\n', ...
        r.stageSeconds, d.nSteps, d.nChannelsAffected, r.nbchan, d.nSteps / r.minutes, ...
        numel(d.events.id), d.medianRiseMs, median(abs(S.steps.amplitude)));
    fprintf('  rejected %d bad contact(s) (%s; before eye flagging the detector named %d: %s; eye electrodes protected: %s); %d common-mode steps in %d events; %d eye-movement steps in %d events (mask %.1f%%: blink %.1f%%, HEOG %.1f%%, %s)\n', ...
        numel(r.rejected), strjoin(r.rejected, ' '), numel(d.badChannelsBeforeEyeFlag), strjoin(d.badChannelsBeforeEyeFlag, ' '), ...
        strjoin(d.protectedEyeChannels', ' '), ...
        d.nCommonModeSteps, d.nCommonModeEvents, ...
        e.nEyeSteps, e.nEyeEvents, 100*e.mask.maskFraction, 100*e.mask.blinkFraction, 100*e.mask.heogFraction, e.mask.heogStatus);
    fprintf('  rejections by detector: amplitude %d, edge %d\n', ...
        d.nRejected.amplitude, d.nRejected.edge);
    fprintf('  repair: %d intervals filled (%d spatial, %d linear, median basis %d), median %.0f ms; %.3f%% of samples, worst channel %.2f%%; %d excised (%.2f s), %d too long; %d sub-segments, %d fallbacks\n', ...
        p.nIntervals, p.nFillsSpatial, p.nFillsLinear, p.medianBasisSize, p.medianTransitionMs, ...
        100*p.fractionSamplesInterpolated, 100*p.maxChannelFractionInterpolated, p.nExcised, p.excisedSeconds, ...
        p.nCommonModeTooLong, p.nSubSegments, p.treatmentCounts.fallback);
    over = find(p.fractionInterpolatedPerChannel > 0.02);
    if ~isempty(over)
        fprintf('  ** %d channel(s) over the 2%% interpolation ceiling: %s\n', numel(over), ...
            strjoin(arrayfun(@(c) sprintf('%s(%.1f%%)', S.preRepairEEG.chanlocs(c).labels, ...
            100*p.fractionInterpolatedPerChannel(c)), over(:)', 'UniformOutput', false), ' '));
    end
end

% ------------------------------------------------------------------------------------------
function R = ringingAtSites(V, F, srate)
% std of the 1 Hz-highpassed channel over [lo - 1 s, hi + 1 s] per fill and version, plus a
% step-free control window on the same channel (same length, half a file away).
    names = fieldnames(V);
    pad = round(srate);
    R = struct('n', 0, 'versions', {names'});
    HP = struct();
    for v = 1:numel(names)
        E = pop_eegfiltnew(V.(names{v}), 1, []);
        HP.(names{v}) = double(E.data);
    end
    n = size(F, 1); pnts = size(HP.raw, 2);
    vals = NaN(n, numel(names)); ctrl = NaN(n, 1); len = F(:,3) - F(:,2) + 1;
    for i = 1:n
        c = F(i,1); a = F(i,2) - pad; b = F(i,3) + pad;
        if a < 1 || b > pnts; continue; end
        for v = 1:numel(names); vals(i, v) = std(HP.(names{v})(c, a:b)); end
        w = b - a;
        j = mod(a + round(pnts/2), pnts - w - 1) + 1;
        ctrl(i) = std(HP.spatial(c, j:j+w));
    end
    ok = ~isnan(ctrl);
    R.n = sum(ok); R.lenMs = len(ok) / srate * 1000;
    R.perFill = vals(ok, :); R.control = ctrl(ok);
    for v = 1:numel(names)
        R.median.(names{v}) = median(vals(ok, v));
        R.ratioToControl.(names{v}) = median(vals(ok, v) ./ ctrl(ok));
    end
    R.controlMedian = median(ctrl(ok));
    long = R.lenMs > 50;
    R.fracSpatialLeqLinear = mean(R.perFill(:, strcmp(names, 'spatial')) <= R.perFill(:, strcmp(names, 'linear')));
    R.fracSpatialLeqLinearLong = mean(R.perFill(long, strcmp(names, 'spatial')) <= R.perFill(long, strcmp(names, 'linear')));
end

function printRinging(R)
    fprintf('post-1 Hz-HP std over [onset-1 s, offset+1 s], %d fills (control %.2f uV):', R.n, R.controlMedian);
    for v = R.versions
        fprintf(' %s %.2f uV (%.2fx)', v{1}, R.median.(v{1}), R.ratioToControl.(v{1}));
    end
    fprintf(['\n  spatial <= linear on %.0f%% of fills (%.0f%% of fills > 50 ms) -- a content comparison, not an ' ...
             'error one: the spatial fill puts neighbour detail back into the window; see the LOO line for error\n'], ...
        100*R.fracSpatialLeqLinear, 100*R.fracSpatialLeqLinearLong);
end

% ------------------------------------------------------------------------------------------
function L = spatialLeaveOneOut(EEG, F, basisExclude, anchorMs)
% For each fill (c, lo, hi) take a step-free window of the same length on the same channel
% and ask how much of the channel's actual anchor-relative detail the neighbours' spline
% estimate explains. Same basis rule as RELAX_repair_steps (coordinates, not c, not
% excluded, not busy), same anchors, same weights.
    L = struct('n', 0, 'rmsLinear', NaN, 'rmsSpatial', NaN, 'fracSpatialBetter', NaN, 'fracSpatialBetterLong', NaN);
    Y = double(EEG.data); [nb, pnts] = size(Y);
    k = max(1, round(anchorMs * EEG.srate / 1000));
    hasXYZ = false(nb, 1); XYZ = NaN(nb, 3);
    for c = 1:nb
        if all(isfield(EEG.chanlocs, {'X','Y','Z'})) && ~isempty(EEG.chanlocs(c).X)
            XYZ(c,:) = [EEG.chanlocs(c).X EEG.chanlocs(c).Y EEG.chanlocs(c).Z]; hasXYZ(c) = all(isfinite(XYZ(c,:)));
        end
    end
    excluded = ismember({EEG.chanlocs.labels}, basisExclude)';
    segments = RELAX_get_segments(EEG);
    rmsLin = []; rmsSp = []; lens = [];
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:size(F, 1)
        c = F(i,1); len = F(i,3) - F(i,2) + 1;
        if ~hasXYZ(c); continue; end
        % control window: shift by half the file, then move to a spot inside one segment
        % that no fill on any channel touches (so the basis is the full clean montage)
        j = mod(F(i,2) + round(pnts/2), pnts - len - 2*k - 1) + k + 1;
        found = false;
        for attempt = 1:20
            lo = j; hi = j + len - 1;
            inSeg = any(segments(:,1) <= lo - k & segments(:,2) >= hi + k);
            busy = F(:,3) >= lo - k - 1 & F(:,2) <= hi + k + 1;
            if inSeg && ~any(busy); found = true; break; end
            j = mod(j + round(pnts/37), pnts - len - 2*k - 1) + k + 1;
        end
        if ~found; continue; end
        B = find(hasXYZ & ~excluded); B(B == c) = [];
        if numel(B) < 16; continue; end
        key = sprintf('%d:%s', c, mat2str(B'));
        if isKey(cache, key); w = cache(key); else; w = RELAX_spline_weights(XYZ(B,:), XYZ(c,:)); cache(key) = w; end
        rows = [c; B];
        Lm = mean(Y(rows, lo-k:lo-1), 2); Rm = mean(Y(rows, hi+1:hi+k), 2);
        chord = Lm + (Rm - Lm) .* ((1:len) / (len + 1));
        detail = Y(rows, lo:hi) - chord;
        est = w * detail(2:end, :);
        rmsLin(end+1) = rms(detail(1,:)); %#ok<AGROW>
        rmsSp(end+1)  = rms(detail(1,:) - est); %#ok<AGROW>
        lens(end+1)   = len / EEG.srate * 1000; %#ok<AGROW>
    end
    L.n = numel(rmsLin);
    if L.n == 0; return; end
    L.rmsLinear = median(rmsLin); L.rmsSpatial = median(rmsSp);
    L.fracSpatialBetter = mean(rmsSp <= rmsLin);
    L.fracSpatialBetterLong = mean(rmsSp(lens > 50) <= rmsLin(lens > 50));
    L.perWindow = [lens(:) rmsLin(:) rmsSp(:)];
end

% ------------------------------------------------------------------------------------------
function r = compareICA(V, nPCA)
% PICARD on every version in V (V.raw is the reference) after the pipeline's own
% preparation, on the channels all versions share; ICs of each other version matched
% one-to-one to raw's by |corr| of the mixing columns, reported overall and per ICLabel
% class of the raw IC. Also the floor of the statistic: raw against raw minus its last 1 %
% of samples, i.e. how much the decomposition moves under a trivial change of the data.
    names = fieldnames(V);
    labels = {V.raw.chanlocs.labels};
    for k = 1:numel(names); labels = intersect(labels, {V.(names{k}).chanlocs.labels}, 'stable'); end
    prep = @(E) pop_resample(pop_eegfiltnew(pop_select(E, 'channel', labels), 1, []), 250);
    ica  = @(E) pop_iclabel(pop_runica(E, 'icatype', 'picard', 'pca', nPCA, 'maxiter', 500), 'default');
    classes = {'Brain','Muscle','Eye','Heart','LineNoise','ChanNoise','Other'};
    P = struct(); secs = struct();
    for k = 1:numel(names)
        E = prep(V.(names{k})); tic; P.(names{k}) = ica(E); secs.(names{k}) = toc;
    end
    E0 = prep(V.raw); E0 = pop_select(E0, 'point', [1 round(0.99 * E0.pnts)]);
    P.floor = ica(E0); names{end+1} = 'floor';
    [~, classRaw] = max(P.raw.etc.ic_classification.ICLabel.classifications, [], 2);
    r = struct('nChannels', numel(labels), 'nPCA', nPCA, 'classNames', {classes}, ...
               'classRaw', classRaw, 'classHistRaw', accumarray(classRaw, 1, [7 1])', 'seconds', secs);
    for k = 2:numel(names)
        nm = names{k};
        [rr, M] = matchComponents(P.raw.icawinv, P.(nm).icawinv);
        [~, cls] = max(P.(nm).etc.ic_classification.ICLabel.classifications, [], 2);
        c = struct('nPairs', size(M,1), 'medianAbsR', median(rr), 'nAbove090', sum(rr > 0.9), ...
                   'nAbove095', sum(rr > 0.95), 'matchedAbsR', rr, 'pairs', M, 'class', cls, ...
                   'classHist', accumarray(cls, 1, [7 1])');
        % matched |r| by the class ICLabel gave the RAW component of each pair
        c.medianAbsRByRawClass = NaN(1, 7); c.nByRawClass = zeros(1, 7);
        for j = 1:7
            sel = classRaw(M(:,1)) == j; c.nByRawClass(j) = sum(sel);
            if any(sel); c.medianAbsRByRawClass(j) = median(rr(sel)); end
        end
        r.(nm) = c;
    end
    % every pair of arms, not only each against raw
    r.arms = names'; nA = numel(names);
    r.pairMedianAbsR = NaN(nA); r.pairNAbove095 = NaN(nA);
    for i = 1:nA
        for j = 1:nA
            rr = matchComponents(P.(names{i}).icawinv, P.(names{j}).icawinv);
            r.pairMedianAbsR(i, j) = median(rr); r.pairNAbove095(i, j) = sum(rr > 0.95);
        end
    end
    % split-half reproducibility: a reference that needs no pop-free data. If the repair
    % removes non-stationary artifact structure, the two halves agree more after it.
    r.splitHalf = struct();
    for nm = {'raw', 'repaired'}
        E = prep(V.(nm{1})); h = floor(E.pnts / 2);
        H1 = pop_runica(pop_select(E, 'point', [1 h]), 'icatype', 'picard', 'pca', nPCA, 'maxiter', 500);
        H2 = pop_runica(pop_select(E, 'point', [h+1 E.pnts]), 'icatype', 'picard', 'pca', nPCA, 'maxiter', 500);
        rr = matchComponents(H1.icawinv, H2.icawinv);
        r.splitHalf.(nm{1}) = struct('medianAbsR', median(rr), 'nAbove090', sum(rr > 0.9), ...
                                     'nAbove095', sum(rr > 0.95), 'matchedAbsR', rr);
    end
    r.classSummary = strjoin(arrayfun(@(k) sprintf('%s %d/%d', classes{k}, r.classHistRaw(k), ...
        r.repaired.classHist(k)), 1:7, 'UniformOutput', false), ', ');
    r.byClassSummary = strjoin(arrayfun(@(k) sprintf('%s %.2f (n=%d)', classes{k}, ...
        r.repaired.medianAbsRByRawClass(k), r.repaired.nByRawClass(k)), 1:7, 'UniformOutput', false), ', ');
end

function [rr, M] = matchComponents(W1, W2)
% One-to-one matching of mixing columns by |corr| (Hungarian assignment on 1 - |corr|).
    C = abs(corr(W1, W2));
    M = matchpairs(1 - C, 1);                         % rows [i1 i2]
    rr = C(sub2ind(size(C), M(:,1), M(:,2)));
end

% ------------------------------------------------------------------------------------------
function Lg = legacyTools(enable)
% The previous detector/repair pair (commit f6649a8) under *_legacy names, for the
% same-site comparison. Empty if git or the checkout is unavailable. --no-pager is
% required: MATLAB's system() gives the child a pseudo-terminal, so git would otherwise
% hand the blob to less and block forever.
    Lg = [];
    if ~enable; return; end
    repo = fileparts(mfilename('fullpath'));
    out = fullfile(tempdir, 'relax_pop_legacy');
    if ~exist(out, 'dir'); mkdir(out); end
    for name = {'RELAX_detect_steps', 'RELAX_repair_steps'}
        [st, src] = system(sprintf('git --no-pager -C "%s" show f6649a8:%s.m', repo, name{1}));
        if st ~= 0
            warning('legacy %s unavailable (git show failed); skipping the legacy comparison', name{1});
            return
        end
        src = strrep(src, name{1}, [name{1} '_legacy']);
        fid = fopen(fullfile(out, [name{1} '_legacy.m']), 'w'); fwrite(fid, src); fclose(fid);
    end
    addpath(out);
    Lg = struct('detect', @(E) RELAX_detect_steps_legacy(E), 'repair', @(E, s) RELAX_repair_steps_legacy(E, s));
end

% ------------------------------------------------------------------------------------------
function printSummary(results)
    fprintf('\n===== SUMMARY over %d files\n', numel(results));
    fprintf('%-12s %6s %6s %5s %4s %5s %6s %6s %7s %6s %6s %6s %6s %6s\n', 'file', 'steps', '/min', 'rej', 'CM', 'eye', 'fills', 'exc_s', 'maxCh%', 'raw', 'legacy', 'linear', 'spat', 'ctrl');
    for i = 1:numel(results)
        r = results{i}; rg = r.ringing; leg = NaN; if isfield(rg.median, 'legacy'); leg = rg.median.legacy; end
        fprintf('%-12s %6d %6.1f %5d %4d %5d %6d %6.2f %7.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n', ...
            r.file(1:min(12,end)), r.detect.nSteps, r.detect.nSteps / r.minutes, numel(r.rejected), ...
            r.detect.nCommonModeEvents, r.eye.nEyeSteps, r.repair.nIntervals, r.repair.excisedSeconds, ...
            100*r.repair.maxChannelFractionInterpolated, rg.median.raw, leg, rg.median.linear, ...
            rg.median.spatial, rg.controlMedian);
    end
    fprintf('(raw/legacy/linear/spatial/ctrl: median std in uV of the 1 Hz-highpassed channel around each filled transition)\n');
end

function s = setDefault(s, f, v)
    if ~isfield(s, f); s.(f) = v; end
end
