% RELAX_bridged_channels() - Detect electrically bridged electrode clusters via eBridge.
%
% Usage:
%   >> [rejectLabels, info] = RELAX_bridged_channels(EEG);
%   >> [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1, 'BCT', 0.5);
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .chanlocs, .srate, .nbchan, .pnts, .event)
%   Name-value pairs (all optional):
%     EpochSec              - epoch length in seconds for the electrical distance
%                              computation (default 1; eBridge's own auto-range is
%                              0.5-1.0 s)
%     Filter                 - 1 = bandpass a working copy before computing ED, via
%                              RELAX_filtbutter; 0 = use the data as-is (default 0).
%                              Off by default: band limiting rescales every ED value,
%                              and this screen normally runs on data already notched.
%     FiltLow, FiltHigh       - bandpass edges in Hz, used only if Filter=1 (default 0.5, 30)
%     FiltOrder               - Butterworth order for the pre-ED bandpass (default 2)
%     BCT                     - eBridge BCT: fraction of epochs below the ED cutoff
%                              required to call a pair bridged (default 0.5)
%     BinSize                 - eBridge BinSize: bin width of the ED frequency
%                              distribution (default 0.1)
%     MaxEpochs                - cap on epochs handed to eBridge, matching eBridge's own
%                              internal cap (default 1500); subsampled deterministically
%     MaxProportionChannels    - QC only, does not change what is rejected: proportion of
%                              channels above which bridging is flagged as excessive
%                              (default 0.10)
%     MaxClusterSize          - QC only: cluster size above which interpolation over the
%                              resulting gap is considered untrustworthy (default 4)
%     RespectBoundaries        - passed through to RELAX_filtbutter when Filter=1 (default 0)
%
% Outputs:
%   rejectLabels - Column cell array of channel labels to reject: every channel
%                  belonging to a bridged cluster (all members, not all-but-one - a
%                  bridge is a physical short and the two original signals are not
%                  recoverable from one measurement; a retained member would still hold
%                  the smeared mixture while its chanlocs entry claims a precise
%                  location). {} if none / eBridge unavailable / too few epochs.
%   info          - struct of QC fields (store on EEG.RELAX for cohort-level review):
%                     .available          false when eBridge is missing, there were too
%                                         few usable epochs, or eBridge itself errored;
%                                         nothing is rejected
%                     .status             human-readable reason when not available
%                     .filtered, .filtBand, .epochSec, .numEpochs
%                     .EDcutoff, .EDscale
%                     .clusters, .pairs, .nBridgedChannels, .nBridgedPairs
%                     .oversizedCluster, .excessiveBridging   QC flags, see above -
%                                         these do not change what is rejected
%                     .rejectLabels
%
% Computes the pairwise electrical distance ED(i,j) = var(x_i - x_j) over short epochs
% and flags pairs whose ED sits in the near-zero peak of the ED distribution (Alschuler
% et al. 2014, Clin Neurophysiol 125:484-490), via the EEGLAB ebridge plugin's eBridge()
% (not vendored in this repo - install the ebridge plugin to use this screen). ED is
% reference-invariant, so it is valid regardless of the current reference. It is also
% amplitude-sensitive, which matters on high-density nets where adjacent electrodes
% routinely correlate at r > 0.95 for purely neurophysiological reasons - a correlation
% criterion would flag them, ED does not.
%
% Epochs are built from RELAX_get_segments so they respect boundary events instead of
% straddling excised (e.g. TMS) artifact windows; a plain reshape would silently corrupt
% ED on any file with boundaries in it. eBridge's own FiltMode is not used regardless (it
% is only honoured for continuous input and silently ignored for epoched input), so
% filtering, if wanted, happens here via RELAX_filtbutter before epoching.
%
% Unlike RELAX_impedance_bad_channels, there is no enforced cap on how many channels this
% screen can reject: bridging is measured directly from the data (a robust statistical
% test already ran), not metadata that could itself be wrong, so there is no argument for
% under-removing a channel this test flags. MaxProportionChannels/MaxClusterSize are
% QC-only: they mark a recording whose reconstruction should not be trusted (nothing
% derived from a bridged channel is fed back in, so a large contiguous cluster leaves a
% gap with no interior information, and spherical interpolation will happily produce a
% smooth, confident, meaningless surface across it) without changing what gets rejected.
%
% Ported from detect_bridged_channels / remove_bridged_channels in
% TMS_response_variability/iTBS_pipeline.m; grouping of flagged pairs into clusters is
% split out to RELAX_group_bridged_channels for independent unit testing.

function [rejectLabels, info] = RELAX_bridged_channels(EEG, varargin)

    p = inputParser();
    p.KeepUnmatched = true;
    addParameter(p, 'EpochSec', 1);
    addParameter(p, 'Filter', 0);
    addParameter(p, 'FiltLow', 0.5);
    addParameter(p, 'FiltHigh', 30);
    addParameter(p, 'FiltOrder', 2);
    addParameter(p, 'BCT', 0.5);
    addParameter(p, 'BinSize', 0.1);
    addParameter(p, 'MaxEpochs', 1500);
    addParameter(p, 'MaxProportionChannels', 0.10);
    addParameter(p, 'MaxClusterSize', 4);
    addParameter(p, 'RespectBoundaries', 0);
    parse(p, varargin{:});
    epochSec          = p.Results.EpochSec;
    doFilter          = p.Results.Filter;
    filtLow           = p.Results.FiltLow;
    filtHigh          = p.Results.FiltHigh;
    filtOrder         = p.Results.FiltOrder;
    bct               = p.Results.BCT;
    binSize           = p.Results.BinSize;
    maxEpochs         = p.Results.MaxEpochs;
    maxPropChans      = p.Results.MaxProportionChannels;
    maxClusterSz      = p.Results.MaxClusterSize;
    respectBoundaries = p.Results.RespectBoundaries;

    minEpochs = 30; % eBridge asserts on this internally; check first for a clearer message

    rejectLabels = {};
    info = struct('available', false, 'status', 'ok', ...
                  'filtered', logical(doFilter), 'filtBand', [filtLow, filtHigh], ...
                  'epochSec', epochSec, 'numEpochs', 0, ...
                  'EDcutoff', NaN, 'EDscale', NaN, ...
                  'clusters', {{}}, 'pairs', {{}}, ...
                  'nBridgedChannels', 0, 'nBridgedPairs', 0, ...
                  'oversizedCluster', false, 'excessiveBridging', false, ...
                  'rejectLabels', {{}});

    if isempty(which('eBridge'))
        warning('RELAX_bridged_channels:eBridgeNotFound', ...
            ['eBridge.m not found on the MATLAB path. Skipping bridged channel ' ...
             'detection. Install the EEGLAB ebridge plugin to use this screen.']);
        info.status = 'eBridge not found';
        return;
    end

    if doFilter
        workingEEG = RELAX_filtbutter(EEG, filtLow, filtHigh, filtOrder, 'bandpass', 'acausal', respectBoundaries);
    else
        workingEEG = EEG;
    end

    segments = RELAX_get_segments(workingEEG);
    windows = local_epoch_windows(segments, workingEEG.srate, epochSec);
    numEpochs = size(windows, 1);

    if numEpochs < minEpochs
        warning('RELAX_bridged_channels:tooFewEpochs', ...
            ['Only %d epoch(s) of length %.2f s available (eBridge requires at ' ...
             'least %d). Skipping bridged channel detection.'], numEpochs, epochSec, minEpochs);
        info.status = sprintf('too few epochs (%d < %d)', numEpochs, minEpochs);
        return;
    end

    if numEpochs > maxEpochs
        % Deterministic even subsampling (not random) to keep runs reproducible:
        keep = unique(round(linspace(1, numEpochs, maxEpochs)));
        windows = windows(keep, :);
        numEpochs = size(windows, 1);
    end
    info.numEpochs = numEpochs;

    numPts = windows(1,2) - windows(1,1) + 1;
    epochData = zeros(workingEEG.nbchan, numPts, numEpochs);
    for k = 1:numEpochs
        epochData(:,:,k) = double(workingEEG.data(:, windows(k,1):windows(k,2)));
    end

    % eBridge only reads .data, .chanlocs and .srate, so a minimal struct suffices:
    EPO = struct('data', epochData, 'chanlocs', EEG.chanlocs, 'srate', EEG.srate);
    try
        EB = eBridge(EPO, 'BCT', bct, 'BinSize', binSize, 'PlotMode', 0, 'Verbose', 1);
    catch ebException
        warning('RELAX_bridged_channels:eBridgeFailed', ...
            'eBridge failed (%s). Skipping bridged channel detection.', ebException.message);
        info.status = sprintf('eBridge error: %s', ebException.message);
        return;
    end
    clear epochData EPO

    info.available = true;
    if isfield(EB, 'Info')
        if isfield(EB.Info, 'EDcutoff'), info.EDcutoff = EB.Info.EDcutoff; end
        if isfield(EB.Info, 'EDscale'),  info.EDscale  = EB.Info.EDscale;  end
    end

    % EDcutoff == -1 means no near-zero peak was found in the ED distribution, i.e. no
    % bridging - distinct from an error:
    if info.EDcutoff == -1
        info.status = 'no ED cutoff (no bridging)';
        return;
    end
    if ~isfield(EB, 'Bridged') || EB.Bridged.Count == 0
        return;
    end

    labels = {EEG.chanlocs.labels};
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, EB.Bridged.Pairs);
    if isempty(clusters)
        return;
    end

    clusterSizes = cellfun(@numel, clusters);
    info.clusters = clusters;
    info.pairs = pairLabels;
    info.nBridgedPairs = numel(pairLabels);
    info.nBridgedChannels = sum(clusterSizes);

    memberIdx = [];
    for c = 1:numel(clusters)
        memberIdx = [memberIdx, find(ismember(labels, clusters{c}))]; %#ok<AGROW>
    end
    memberIdx = unique(memberIdx);
    rejectLabels = labels(memberIdx)';
    info.rejectLabels = rejectLabels;

    if any(clusterSizes > maxClusterSz)
        info.oversizedCluster = true;
        warning('RELAX_bridged_channels:oversizedCluster', ...
            ['%d bridged cluster(s) exceed MaxClusterSize (%d); largest is %d ' ...
             'channels. Interpolation over a patch this large is unreliable - ' ...
             'consider excluding this recording.'], ...
            sum(clusterSizes > maxClusterSz), maxClusterSz, max(clusterSizes));
    end
    if info.nBridgedChannels > maxPropChans * numel(labels)
        info.excessiveBridging = true;
        warning('RELAX_bridged_channels:excessiveBridging', ...
            ['%d of %d channels (%.1f%%) are bridged, above the %.0f%% threshold. ' ...
             'Consider excluding this recording.'], ...
            info.nBridgedChannels, numel(labels), 100*info.nBridgedChannels/numel(labels), 100*maxPropChans);
    end
end

function windows = local_epoch_windows(segments, fs, epochSec)
    % Non-overlapping windows of length epochSec within each boundary-respecting
    % segment, never crossing a segment edge.
    windowSamples = round(epochSec * fs);
    windows = zeros(0, 2);
    for s = 1:size(segments,1)
        segStart = segments(s,1);
        segEnd = segments(s,2);
        segLen = segEnd - segStart + 1;
        if segLen < windowSamples
            continue; % segment too short to fit even one window
        end
        winStarts = (segStart:windowSamples:(segEnd - windowSamples + 1))';
        if isempty(winStarts)
            continue;
        end
        winEnds = winStarts + windowSamples - 1;
        windows = [windows; winStarts, winEnds]; %#ok<AGROW>
    end
end
