function test_RELAX_blinks_IQR_fallback()
% RELAX_blinks_IQR_method must not error when fewer than two of RELAX_cfg.BlinkElectrodes
% survive earlier channel rejection. The native extreme-value screen legitimately removed
% 6 of the 8 frontal electrodes on m160 day 2 (they were extreme in too many epochs) and
% the bridge screen removed all 8 on m396; both files then died with "Insufficient blink
% electrodes available". Expected behaviour: fall back to the row of electrodes one row
% higher (E12 E20 E27 E34 E21), and if that still leaves fewer than two, report no blinks
% through the no-blink path every later stage already handles. Either way the choice is
% recorded as a per-file covariate, because a file cleaned without blink marking is
% cleaned differently and the downstream model has to be able to see that.

    if exist('pop_select', 'file') ~= 2
        fprintf('test_RELAX_blinks_IQR_fallback: SKIPPED (EEGLAB not on the path)\n');
        return
    end

    primary = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
    backup  = {'E12'; 'E20'; 'E27'; 'E34'; 'E21'};
    others  = {'E70'; 'E80'; 'E90'; 'E100'; 'E110'; 'E120'; 'E130'; 'E140'; 'E150'; 'E160'};
    srate = 1000;

    cfg = struct();
    cfg.BlinkElectrodes = primary;
    cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks = 'no';
    cfg.ms_per_sample = 1000 / srate;
    cfg.HighPassFilter = 1; cfg.LowPassFilter = 200;
    cfg.RespectBoundaries = 1;

    % 1. every primary electrode present: unchanged behaviour, covariates say so
    EEG = makeEEG([primary; backup; others], srate);
    [c, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(c.RELAX.IQRmethodDetectedBlinks == 1, 'blink is detected on the full montage');
    assert(c.RELAX.eyeblinkmask(10*srate) == 1, 'mask covers the blink at 10 s');
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 8, ...
        'all 8 primaries used, got %d', e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 0, 'no fallback');
    assert(c.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 0, ...
        'covariates are on the continuous output too');

    % 2. one primary left (m160 day 2): backup row joins it, blink still found
    lastwarn('');
    EEG = makeEEG([primary(1); backup; others], srate);
    [c, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(c.RELAX.IQRmethodDetectedBlinks == 1, 'blink is detected from the backup row');
    assert(c.RELAX.eyeblinkmask(10*srate) == 1, 'mask covers the blink at 10 s (backup row)');
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 6, ...
        '1 primary + 5 backups used, got %d', e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 1, 'fallback = backup row');
    [~, wid] = lastwarn();
    assert(strcmp(wid, 'RELAX_blinks_IQR_method:backupElectrodes'), 'a warning names the fallback');

    % 3. no primaries at all (m396): backup row alone
    EEG = makeEEG([backup; others], srate);
    [c, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(c.RELAX.IQRmethodDetectedBlinks == 1, 'blink is detected from backups alone');
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 5);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 1);

    % 4. fewer than two even with the backups: no error, no blinks, all-zero mask
    lastwarn('');
    EEG = makeEEG([primary(1); others], srate);
    [c, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(c.RELAX.IQRmethodDetectedBlinks == 0, 'no blinks reported');
    assert(e.RELAX.IQRmethodDetectedBlinks == 0, 'epoched copy agrees');
    assert(isequal(size(c.RELAX.eyeblinkmask), [1 EEG.pnts]) && ~any(c.RELAX.eyeblinkmask), ...
        'mask is 1 x pnts and all zero');
    assert(isequal(e.RELAX.eyeblinkmask, c.RELAX.eyeblinkmask), 'epoched mask matches');
    assert(isempty(c.event) || ~any(strcmp({c.event.type}, 'EyeBlinkMax')), 'no blink events written');
    assert(isequal(c.data, EEG.data), 'data untouched');
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 0);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 2, 'fallback = none');
    [~, wid] = lastwarn();
    assert(strcmp(wid, 'RELAX_blinks_IQR_method:insufficientElectrodes'), 'a warning names the skip');

    % 5. the 2-electrode threshold is inclusive: two primaries need no fallback, and one
    %    primary plus one backup is enough to run on the backup path
    EEG = makeEEG([primary(1:2); others], srate);
    [~, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 2);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 0);
    EEG = makeEEG([primary(1); backup(1); others], srate);
    [c, e] = evalc_quiet(@() RELAX_blinks_IQR_method(EEG, EEG, cfg));
    assert(c.RELAX.IQRmethodDetectedBlinks == 1, 'one primary + one backup still detects');
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodesUsed == 2);
    assert(e.RELAXProcessingExtremeRejections.BlinkElectrodeFallback == 1);

    fprintf('test_RELAX_blinks_IQR_fallback: PASS\n');
end

function EEG = makeEEG(labels, srate)
% 30 s, one 150 uV blink at 10 s on every E1x-E3x label (the whole periocular region),
% 20 uV elsewhere, over pink-ish background noise. Average re-reference is marked as
% already done so RELAX_average_rereference leaves the data alone.
    pnts = 30 * srate;
    nb = numel(labels);
    rs = RandStream('twister', 'Seed', 3);
    t = (0:pnts-1) / srate;
    data = 5 * randn(rs, nb, pnts) + 20 * sin(2*pi*0.2*t + 2*pi*rand(rs, nb, 1));
    blink = exp(-0.5 * ((t - 10) / 0.06).^2);
    num = cellfun(@(s) sscanf(s, 'E%d'), labels);
    periocular = num <= 40;
    data(periocular, :)  = data(periocular, :)  + 150 * blink;
    data(~periocular, :) = data(~periocular, :) + 20 * blink;

    EEG = eeg_emptyset();
    EEG.setname = 'blinkfallback'; EEG.srate = srate; EEG.nbchan = nb; EEG.pnts = pnts;
    EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (pnts-1)/srate; EEG.data = data;
    xyz = randn(rs, nb, 3); xyz = xyz ./ sqrt(sum(xyz.^2, 2));
    EEG.chanlocs = struct('labels', labels', 'X', num2cell(xyz(:,1))', ...
                          'Y', num2cell(xyz(:,2))', 'Z', num2cell(xyz(:,3))');
    EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
    EEG = eeg_checkset(EEG);
    EEG.allchan = EEG.chanlocs;
    EEG.RELAX = struct('ListOfChannelsAfterRejections', {labels'}, ...
                       'Data_has_been_averagerereferenced', 1);
    EEG.RELAXProcessing = struct();
    EEG.RELAXProcessingExtremeRejections = struct();
end

function varargout = evalc_quiet(fn) %#ok<INUSD>
% Run fn with its command-window chatter captured.
    varargout = cell(1, max(1, nargout));
    [~, varargout{:}] = evalc('fn()');
end
