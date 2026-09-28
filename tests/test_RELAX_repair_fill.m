function test_RELAX_repair_fill()
% Interval repair: transitions are bridged (spatially when neighbours exist), skipped steps
% are left alone, common-mode events are left alone by default and excised with a boundary
% when ExciseCommonMode is on, legacy single-sample
% steps still reproduce the old partition.

    RELAX_add_vendored_noisetools();
    srate = 1000; n = 10000; rng(19);

    function EEG = mkEEG(data, boundaries, XYZ)
        if nargin < 2; boundaries = []; end
        if exist('eeg_emptyset', 'file') == 2
            EEG = eeg_emptyset;                 % full EEGLAB skeleton, so eeg_eegrej can run
        else
            EEG = struct();
        end
        EEG.data = data; EEG.nbchan = size(data,1); EEG.pnts = size(data,2); EEG.srate = srate;
        EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (size(data,2)-1)/srate;
        EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
        for b = boundaries(:)'
            EEG.event(end+1) = struct('type', 'boundary', 'latency', b, 'duration', NaN);
        end
        for c = 1:size(data,1)
            EEG.chanlocs(c).labels = sprintf('E%d', c);
            if nargin >= 3
                EEG.chanlocs(c).X = XYZ(c,1); EEG.chanlocs(c).Y = XYZ(c,2); EEG.chanlocs(c).Z = XYZ(c,3);
            end
        end
    end
    function s = mkSteps(ch, onset, offset, varargin)
        % varargin: 'commonMode', logical vector | 'eyeMovement', logical vector | 'eventId'
        m = numel(onset);
        s = struct('channel', ch(:), 'label', {arrayfun(@(k) sprintf('E%d', k), ch(:), 'UniformOutput', false)}, ...
                   'sample', onset(:) + 1, 'onset', onset(:), 'offset', offset(:), ...
                   'amplitude', 100 * ones(m,1), 'z', 10 * ones(m,1), 'eventId', (1:m)', ...
                   'commonMode', false(m,1), 'eyeMovement', false(m,1));
        for k = 1:2:numel(varargin); s.(varargin{k}) = logical(varargin{k+1}(:)) | false(m,1); end
        if any(strcmp(varargin, 'eventId')); s.eventId = varargin{find(strcmp(varargin,'eventId'))+1}(:); end
    end
    function s = legacySteps(ch, samp)
        s = struct('channel', ch(:), 'label', {cell(numel(ch),1)}, 'sample', samp(:), ...
                   'amplitude', zeros(numel(ch),1), 'z', zeros(numel(ch),1));
    end
    ramp = @(t, a, L) min(1, max(0, (t - a) / L));
    t = 0:n-1;

    %% 1. a ramp transition is bridged linearly when no neighbours exist
    x = randn(1, n) * 5 + 200 * ramp(t, 4000, 150);           % rises over 4001..4150
    [out, info] = RELAX_repair_steps(mkEEG(x), mkSteps(1, 4000, 4151));
    y = out.data(1, :);
    assert(info.available, '%s', info.status);
    assert(abs(median(y(3001:3999)) - median(y(4152:5000))) < 10, 'levels must agree after repair');
    assert(max(abs(diff(y(4000:4151), 2))) < 1e-9, 'the bridge is exactly linear');
    assert(abs(y(4000) - mean(y(3997:3999))) < 1e-9, 'left end sits on the left anchor');
    assert(abs(y(4151) - mean(y(4152:4154))) < 1e-9, 'right end sits on the right anchor');
    assert(~any(isnan(y)));
    assert(info.nIntervals == 1 && isequal(info.intervalsPerChannel, 1));
    assert(abs(info.fractionSamplesInterpolated - 152/n) < 1e-12);
    assert(isequal(info.fractionInterpolatedPerChannel, 152/n));
    assert(info.maxChannelFractionInterpolated == 152/n);
    assert(info.medianTransitionMs == 152);
    assert(info.nFillsLinear == 1 && info.nFillsSpatial == 0, 'no coordinates -> linear bridge');
    assert(info.nCommonModeSkipped == 0 && info.nEyeMovementSkipped == 0 && info.nExcised == 0);

    %% 2. legacy steps (no onset/offset) are zero-length cuts: same partition as before, no fill
    x = randn(1, n) * 5; x(4001:end) = x(4001:end) + 500;
    [a, infoA] = RELAX_repair_steps(mkEEG(x), legacySteps(1, 4001));
    [b, infoB] = RELAX_repair_steps(mkEEG(x), mkSteps(1, 4001, 4000));
    assert(isequal(a.data, b.data), 'legacy sample == explicit [sample, sample-1]');
    assert(infoA.nIntervals == 0 && infoA.legacyIntervals && ~infoB.legacyIntervals);
    assert(abs(median(a.data(1, 4001:5000)) - median(a.data(1, 3001:4000))) < 10);

    %% 3. common-mode (when not excising) and eye-movement steps are left untouched
    x = randn(1, n) * 5 + 200 * ramp(t, 4000, 150);
    ref = RELAX_repair_steps(mkEEG(x), mkSteps(zeros(0,1), zeros(0,1), zeros(0,1)));
    [o3, i3] = RELAX_repair_steps(mkEEG(x), mkSteps(1, 4000, 4151, 'commonMode', true), 'ExciseCommonMode', false);
    assert(isequal(o3.data, ref.data) && i3.nCommonModeSkipped == 1 && i3.subSegmentsPerChannel(1) == 1);
    % not excising is the default: a common-mode step is flagged only - not cut out, not filled,
    % and not a partition point
    [o3, i3] = RELAX_repair_steps(mkEEG(x), mkSteps(1, 4000, 4151, 'commonMode', true));
    assert(isequal(o3.data, ref.data) && o3.pnts == n && i3.nExcised == 0 && isempty(i3.excised) ...
        && i3.nCommonModeSkipped == 1 && i3.nIntervals == 0 && i3.subSegmentsPerChannel(1) == 1, ...
        'by default a common-mode step is left exactly as a run with no steps leaves it');
    [o3, i3] = RELAX_repair_steps(mkEEG(x), mkSteps(1, 4000, 4151, 'eyeMovement', true));
    assert(isequal(o3.data, ref.data) && i3.nEyeMovementSkipped == 1 && i3.nIntervals == 0);

    %% 4. an interval near a segment edge is extended to the edge and held constant
    x = randn(1, n) * 5;
    [o4, i4] = RELAX_repair_steps(mkEEG(x, 6000.5), mkSteps(1, 6005, 6030));
    y = o4.data(1, :);
    assert(i4.subSegmentsPerChannel(1) == 2, 'two segments, the fill touches the edge: got %d', i4.subSegmentsPerChannel(1));
    assert(max(y(6001:6030)) - min(y(6001:6030)) < 1e-9, 'one-sided fill is a constant hold');
    assert(abs(y(6001) - mean(y(6031:6033))) < 1e-9);

    %% 5. an interval covering a whole segment goes through the ladder, not the fill
    x = randn(1, n) * 5 + 300;
    [o5, i5] = RELAX_repair_steps(mkEEG(x, [5000.5 5100.5]), mkSteps(1, 5001, 5100));
    assert(~any(isnan(o5.data(1,:))) && abs(mean(o5.data(1, 5001:5100))) < 5);
    assert(i5.treatmentCounts.meanOnly >= 1 && i5.nIntervals == 0);

    %% 6. overlapping / nearby intervals merge
    x = randn(1, n) * 5;
    [o6, i6] = RELAX_repair_steps(mkEEG(x), mkSteps([1 1], [4000 4020], [4010 4030]));
    assert(i6.nIntervals == 1, 'gap of 9 samples <= MergeSamples merges, got %d', i6.nIntervals);
    assert(max(abs(diff(o6.data(1, 4000:4030), 2))) < 1e-9);

    %% 7. a pulse bracketed by two steps: the plateau becomes a demeaned sub-segment
    x = randn(1, n) * 5; x(4001:4150) = x(4001:4150) + 200;
    [o7, i7] = RELAX_repair_steps(mkEEG(x), mkSteps([1 1], [4000 4150], [4001 4151]));
    assert(abs(median(o7.data(1, 4002:4149)) - median(o7.data(1, 3000:3990))) < 10, 'plateau removed');
    assert(i7.treatmentCounts.meanOnly >= 1);

    %% 8. spatial fill: the missing stretch is rebuilt from the neighbours' detail
    N = 32; k = (0:N-1)'; phi = pi * (3 - sqrt(5)) * k; zc = 1 - k / (N - 1) * 0.9; r = sqrt(1 - zc.^2);
    XYZ = [r .* cos(phi), r .* sin(phi), zc] * 85;
    d = [0.3 0.5 0.8]; d = d / norm(d); gain = (XYZ / 85) * d';          % N x 1 spatial pattern
    src = 30 * sin(2*pi*5*t/srate) + 20 * sin(2*pi*11*t/srate + 1);      % shared source
    X = gain * src + randn(N, n) * 3 + 1000 * randn(N, 1);               % + raw DC offsets
    target = 5;                                                           % interior electrode
    X(target, :) = X(target, :) + 200 * ramp(t, 4000, 150);
    truth = gain(target) * src(4001:4150);
    [oS, iS] = RELAX_repair_steps(mkEEG(X, [], XYZ), mkSteps(target, 4000, 4151));
    [oL, iL] = RELAX_repair_steps(mkEEG(X, [], XYZ), mkSteps(target, 4000, 4151), 'FillMethod', 'linear');
    fillS = oS.data(target, 4001:4150); fillL = oL.data(target, 4001:4150);
    assert(iS.nFillsSpatial == 1 && iS.nFillsLinear == 0 && iS.medianBasisSize == N - 1);
    assert(iL.nFillsSpatial == 0 && iL.nFillsLinear == 1);
    cS = corr(fillS(:), truth(:));
    assert(cS > 0.9, 'spatial fill should track the true signal (r = %.2f)', cS);
    assert(rms(fillS - truth) < 0.5 * rms(fillL - truth), ...
        'spatial residual %.1f should be well below linear residual %.1f', rms(fillS - truth), rms(fillL - truth));
    assert(abs(median(oS.data(target, 3001:3999)) - median(oS.data(target, 4152:5000))) < 10);
    assert(abs(oS.data(target, 4000) - mean(oS.data(target, 3997:3999))) < 3 * 3, 'join stays at anchor-noise level');
    untouched = setdiff(1:N, target);
    assert(isequal(oS.data(untouched, :), oL.data(untouched, :)), 'other channels are unaffected by the fill method');

    %% 9. a co-popping neighbour is excluded from the basis; a listed bad channel too
    [~, i9] = RELAX_repair_steps(mkEEG(X, [], XYZ), mkSteps([target 6], [4000 4010], [4151 4100]));
    assert(i9.nFillsSpatial == 2 && i9.medianBasisSize == N - 2, 'both fills use the other %d channels', N - 2);
    [~, i9b] = RELAX_repair_steps(mkEEG(X, [], XYZ), mkSteps(target, 4000, 4151), 'BasisExcludeLabels', {'E1','E2'});
    assert(i9b.medianBasisSize == N - 3);

    %% 10. too few basis channels -> linear bridge, counted
    [~, i10] = RELAX_repair_steps(mkEEG(X, [], XYZ), mkSteps(target, 4000, 4151), 'MinBasisChannels', 40);
    assert(i10.nFillsLinear == 1 && i10.nFillsSpatial == 0);

    %% 11. with ExciseCommonMode, common-mode events are excised with a boundary event
    % (needs EEGLAB's eeg_eegrej)
    if exist('eeg_eegrej', 'file') == 2
    Y = randn(5, n) * 5; Y(1:4, 5001:end) = Y(1:4, 5001:end) + 300;
    s11 = mkSteps([1 2 3 4], [5000 5000 5000 5000], [5010 5010 5010 5010], 'commonMode', true(4,1), 'eventId', [1 1 1 1]);
    [o11, i11] = RELAX_repair_steps(mkEEG(Y), s11, 'ExciseCommonMode', true);
    assert(o11.pnts == n - 11 && size(o11.data, 2) == n - 11, 'the 11 excised samples are gone');
    bnd = o11.event(strcmp({o11.event.type}, 'boundary'));
    assert(numel(bnd) == 1 && abs(bnd.latency - 4999.5) < 1 && bnd.duration == 11, ...
        'one boundary of duration 11 near 5000 (got latency %g, duration %g)', bnd.latency, bnd.duration);
    assert(isequal(i11.excised, [5000 5010]) && i11.nExcised == 1 && abs(i11.excisedSeconds - 0.011) < 1e-12);
    assert(i11.nCommonModeSkipped == 0 && i11.nIntervals == 0);
    assert(all(i11.subSegmentsPerChannel == 2), 'the event is a hard cut for every channel');
    assert(abs(median(o11.data(1, 5001:5500)) - median(o11.data(1, 4500:4999))) < 10, 'levels agree across the cut');
    % too long to excise -> left alone and counted
    [o11b, i11b] = RELAX_repair_steps(mkEEG(Y), s11, 'ExciseCommonMode', true, 'MaxExciseMs', 5);
    assert(o11b.pnts == n && i11b.nExcised == 0 && i11b.nCommonModeTooLong == 1 && i11b.nCommonModeSkipped == 4);
    % a common-mode event that is also an eye movement is real EOG, not head movement:
    % left in place (not excised, not filled), counted as eye-movement steps
    s11c = mkSteps([1 2 3 4], [5000 5000 5000 5000], [5010 5010 5010 5010], ...
                   'commonMode', true(4,1), 'eyeMovement', true(4,1), 'eventId', [1 1 1 1]);
    [o11c, i11c] = RELAX_repair_steps(mkEEG(Y), s11c, 'ExciseCommonMode', true);
    assert(o11c.pnts == n && i11c.nExcised == 0 && isempty(i11c.excised), 'eye-flagged common-mode event is not excised');
    assert(i11c.nEyeMovementSkipped == 4 && i11c.nCommonModeSkipped == 0 && i11c.nIntervals == 0);
    assert(all(i11c.subSegmentsPerChannel == 1), 'no cut is made for the skipped event');
    % an excision that would leave less than MinSegmentSeconds (1 s) between itself and a
    % segment edge is extended to the edge: a sub-second sliver would be orphaned by the
    % _boundary epoching and a ~1-sample one crashes zapline's pwelch (m292 day 3)
    s11d = mkSteps([1 2 3 4], 9500 * ones(1,4), 9700 * ones(1,4), 'commonMode', true(4,1), 'eventId', [1 1 1 1]);
    [o11d, i11d] = RELAX_repair_steps(mkEEG(Y), s11d, 'ExciseCommonMode', true);
    assert(isequal(i11d.excised, [9500 n]), 'excision within 1 s of the segment end runs to the end (got %s)', mat2str(i11d.excised));
    assert(o11d.pnts == 9499 && abs(i11d.excisedSeconds - 0.501) < 1e-12);
    % two excisions closer than MinSegmentSeconds merge
    s11e = mkSteps([1 2 3 4 1 2 3 4], [3000 3000 3000 3000 3500 3500 3500 3500], [3200 3200 3200 3200 3700 3700 3700 3700], ...
                   'commonMode', true(8,1), 'eventId', [1 1 1 1 2 2 2 2]);
    [o11e, i11e] = RELAX_repair_steps(mkEEG(Y), s11e, 'ExciseCommonMode', true);
    assert(isequal(i11e.excised, [3000 3700]), 'excisions 0.3 s apart merge (got %s)', mat2str(i11e.excised));
    assert(o11e.pnts == n - 701 && i11e.nExcised == 1);
    % MinSegmentSeconds 0 restores the exact cut
    [~, i11f] = RELAX_repair_steps(mkEEG(Y), s11d, 'ExciseCommonMode', true, 'MinSegmentSeconds', 0);
    assert(isequal(i11f.excised, [9500 9700]));
    else
        fprintf('  (eeg_eegrej not on path - skipping the excision check)\n');
    end

    %% 11g. a constant sub-segment is left as it is, never demeaned to exact zero
    % (m544 day 2: a railed channel, constant per segment, came out globally flat and
    % crashed pop_jointprob; the detector now rejects such channels, and the repair must not
    % manufacture a flat channel out of one that is dead in only some segments)
    Z = randn(2, n) * 5; Z(1, 1:5000) = 300; Z(1, 5001:end) = randn(1, 5000) * 5;
    [o11g, i11g] = RELAX_repair_steps(mkEEG(Z, 5000.5), mkSteps(2, 7000, 7050));
    assert(isequal(o11g.data(1, 1:5000), Z(1, 1:5000)), 'constant sub-segment must come back unchanged');
    assert(max(o11g.data(1, :)) - min(o11g.data(1, :)) > 0, 'the channel must not be flat');
    assert(i11g.treatmentCounts.constant == 1 && i11g.treatmentCounts.meanOnly == 0);

    %% 12. nt_detrend missing: every new field present and shaped
    vendorDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'vendor', 'noisetools');
    rmpath(vendorDir);
    cleanupPath = onCleanup(@() RELAX_add_vendored_noisetools());
    [o12, i12] = RELAX_repair_steps(mkEEG(randn(2, n)), mkSteps(1, 4000, 4151));
    assert(~i12.available && isequal(o12.data, o12.data));
    for f = {'nIntervals','fractionSamplesInterpolated','maxChannelFractionInterpolated', ...
             'medianTransitionMs','nCommonModeSkipped','nEyeMovementSkipped','nFillsSpatial', ...
             'nFillsLinear','medianBasisSize','nExcised','excisedSeconds','nCommonModeTooLong','legacyIntervals'}
        assert(isfield(i12, f{1}), 'info.%s missing on the failure path', f{1});
    end
    assert(isequal(size(i12.intervalsPerChannel), [2 1]) && isequal(size(i12.fractionInterpolatedPerChannel), [2 1]));
    assert(isequal(size(i12.excised), [0 2]));
    clear cleanupPath

    fprintf('test_RELAX_repair_fill: PASS\n');
end
