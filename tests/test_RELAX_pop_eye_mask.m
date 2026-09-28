function test_RELAX_pop_eye_mask()
% RELAX_pop_eye_mask builds, on a temporary copy of the raw recording, the eye-movement
% mask the pop stage hands to RELAX_flag_eye_steps: blink periods from
% RELAX_blinks_IQR_method and horizontal eye movements from RELAX_horizontaleye. The raw
% EEG is never modified and every failure degrades to an all-false mask.

    if exist('pop_select', 'file') ~= 2
        fprintf('test_RELAX_pop_eye_mask: SKIPPED (EEGLAB not on the path)\n');
        return
    end

    srate = 1000; pnts = 30 * srate;
    blinkLabels = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
    otherLabels = {'E46'; 'E10'; 'E70'; 'E80'; 'E90'; 'E100'; 'E110'; 'E120'; 'E130'; 'E140'; 'E150'; 'E160'};
    labels = [blinkLabels; otherLabels];
    nb = numel(labels);
    rng(7);
    t = (0:pnts-1) / srate;
    data = 5 * randn(nb, pnts) + 3000 * randn(nb, 1) + 20 * sin(2*pi*0.1*t + 2*pi*rand(nb,1));
    % blink at 10 s: 150 uV Gaussian on the peri-ocular electrodes, 40 uV elsewhere
    blink = exp(-0.5 * ((t - 10) / 0.06).^2);
    data(1:8, :)   = data(1:8, :)   + 150 * blink;
    data(9:end, :) = data(9:end, :) + 40 * blink;
    % saccade at 20 s: +80 uV on E10 (right), -80 uV on E46 (left) for 300 ms
    sac = t >= 20 & t < 20.3;
    data(strcmp(labels, 'E10'), sac) = data(strcmp(labels, 'E10'), sac) + 80;
    data(strcmp(labels, 'E46'), sac) = data(strcmp(labels, 'E46'), sac) - 80;

    EEG = eeg_emptyset();
    EEG.setname = 'eyemask'; EEG.srate = srate; EEG.nbchan = nb; EEG.pnts = pnts;
    EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (pnts-1)/srate; EEG.data = data;
    % coordinates on the unit sphere so the average-reference path (pop_interp) runs
    xyz = randn(nb, 3); xyz = xyz ./ sqrt(sum(xyz.^2, 2));
    EEG.chanlocs = struct('labels', labels', 'X', num2cell(xyz(:,1))', ...
                          'Y', num2cell(xyz(:,2))', 'Z', num2cell(xyz(:,3))');
    EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
    EEG = eeg_checkset(EEG);
    EEG.allchan = EEG.chanlocs;

    cfg = struct();
    cfg.BlinkElectrodes = blinkLabels;
    cfg.HEOGLeftpattern  = ["E46", "E54"];
    cfg.HEOGRightpattern = ["E10", "E1"];
    cfg.RespectBoundaries = 1;
    cfg.HorizontalEyeMovementType = 2;
    cfg.HorizontalEyeMovementThreshold = 2;
    cfg.HorizontalEyeMovementTimepointsExceedingThreshold = 25;
    cfg.HorizontalEyeMovementTimepointsTestWindow = 49;
    cfg.HorizontalEyeMovementFocus = 200;
    cfg.HighPassFilter = 1; cfg.LowPassFilter = 200;
    cfg.PopEye_IncludeHEOG = 1;

    % 1. both detectors fire where they should, nowhere else
    before = EEG;
    [mask, info] = RELAX_pop_eye_mask(EEG, cfg);
    assert(islogical(mask) && isequal(size(mask), [1 pnts]), 'mask is logical 1 x pnts');
    assert(mask(10*srate), 'blink at 10 s is inside the mask');
    assert(mask(20.15*srate), 'saccade at 20 s is inside the mask');
    assert(~mask(5*srate) && ~mask(25*srate) && ~mask(15*srate), 'quiet periods are outside');
    assert(mean(mask) < 0.10, 'mask covers a small fraction of the file, got %.3f', mean(mask));
    assert(strcmp(info.blinkStatus, 'ok') && strcmp(info.heogStatus, 'ok'));
    assert(info.blinkFraction > 0 && info.heogFraction > 0);
    assert(isequal(before, EEG), 'the input EEG is not modified');
    assert(info.rereferenced, 'the temporary copy was average re-referenced');

    % 2. blink component only: the saccade must come from the HEOG detector
    cfgNoHeog = cfg; cfgNoHeog.HEOGLeftpattern = ["Q1"]; cfgNoHeog.HEOGRightpattern = ["Q2"];
    [m2, info2] = RELAX_pop_eye_mask(EEG, cfgNoHeog);
    assert(m2(10*srate) && ~m2(20.15*srate), 'without HEOG channels only the blink is masked');
    assert(~strcmp(info2.heogStatus, 'ok') && strcmp(info2.blinkStatus, 'ok'));
    assert(info2.heogFraction == 0);

    % 3. no blink electrodes at all: no error, all-false blink part, HEOG still works
    cfgNoBlink = cfg; cfgNoBlink.BlinkElectrodes = {'Q1'; 'Q2'};
    [m3, info3] = RELAX_pop_eye_mask(EEG, cfgNoBlink);
    assert(~m3(10*srate) && m3(20.15*srate));
    assert(~strcmp(info3.blinkStatus, 'ok') && info3.blinkFraction == 0);

    % 4. nothing available -> all false, still no error
    cfgNone = cfgNoBlink; cfgNone.HEOGLeftpattern = ["Q1"]; cfgNone.HEOGRightpattern = ["Q2"];
    [m4, info4] = RELAX_pop_eye_mask(EEG, cfgNone);
    assert(~any(m4) && islogical(m4) && isequal(size(m4), [1 pnts]));
    assert(~strcmp(info4.blinkStatus, 'ok') && ~strcmp(info4.heogStatus, 'ok'));

    % 5. the HEOG contribution is off unless PopEye_IncludeHEOG is set: the saccade is not
    %    masked, but the pair is still resolved so the wrapper can count it as eye
    %    electrodes (protected from rejection, needed by RELAX_horizontaleye later)
    cfgDefault = rmfield(cfg, 'PopEye_IncludeHEOG');
    [m5, info5] = RELAX_pop_eye_mask(EEG, cfgDefault);
    assert(m5(10*srate) && ~m5(20.15*srate), 'by default only the blink detector contributes');
    assert(strcmp(info5.heogStatus, 'disabled') && info5.heogFraction == 0);
    assert(isequal(info5.heogChannels, {'E46', 'E10'}), 'the HEOG pair is resolved even when its detector is off');
    assert(strcmp(info5.blinkStatus, 'ok'));

    fprintf('test_RELAX_pop_eye_mask: PASS\n');
end
