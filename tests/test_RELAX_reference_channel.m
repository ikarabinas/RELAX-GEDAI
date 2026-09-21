function test_RELAX_reference_channel()
% Virtual reference channel: the cross-channel median of per-segment-demeaned data.
% It recovers a common (reference) waveform, ignores DC offsets and a minority of channels
% with their own artifacts, is demeaned per inter-boundary segment, and leaves out flat or
% excluded channels. Tolerances: with N = 64 channels of sd-10 noise the median's sd is
% about 1.25 * 10 / 8 = 1.6 uV, and the max over 10^4 samples of that is about 7 uV.

    srate = 1000; n = 10000; rng(3);

    function EEG = mkEEG(data, boundaries)
        if nargin < 2; boundaries = []; end
        EEG = struct();
        EEG.data = data; EEG.nbchan = size(data,1); EEG.pnts = size(data,2);
        EEG.srate = srate; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (EEG.pnts-1)/srate;
        EEG.event = struct('type', {}, 'latency', {});
        for b = boundaries(:)'
            EEG.event(end+1) = struct('type', 'boundary', 'latency', b);
        end
        EEG.chanlocs = struct('labels', arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end
    demean = @(v) v - median(v);

    %% 1. shape and fields
    EEG = mkEEG(randn(8, n), 4000.5);
    [ref, info] = RELAX_reference_channel(EEG);
    assert(info.available, 'available on ordinary data: %s', info.status);
    assert(isequal(size(ref.data), [1 n]) && isa(ref.data, 'double'));
    assert(ref.nbchan == 1 && ref.pnts == n && ref.srate == srate && ref.trials == 1);
    assert(numel(ref.chanlocs) == 1 && strcmp(ref.chanlocs(1).labels, 'REFMEDIAN'));
    assert(numel(ref.event) == 1 && strcmp(ref.event(1).type, 'boundary') && ref.event(1).latency == 4000.5, ...
        'boundary events are carried over so RELAX_detect_steps sees the same segments');
    assert(info.nChannelsUsed == 8 && isequal(size(info.segments), [2 2]));
    assert(iscell(info.usedLabels) && isequal(size(info.usedLabels), [8 1]));
    assert(iscell(info.excludedLabels) && isequal(size(info.excludedLabels), [0 1]));

    %% 2. recovers a common waveform subtracted from every channel
    N = 64; noise = 10;
    r = zeros(1, n); r(5000:end) = 80;                       % the reference pops by +80 uV
    r(4990:5010) = 80 ./ (1 + exp(-(-10:10)/3));             % over ~20 ms
    S = noise * randn(N, n);
    X = S - r;                                                 % x_i = s_i - r
    [ref, info] = RELAX_reference_channel(mkEEG(X));
    err = demean(ref.data) - demean(-r);
    assert(max(abs(err)) < 10, 'median over 64 channels recovers -r to within 10 uV (got %.1f)', max(abs(err)));
    [~, infoS] = RELAX_reference_channel(mkEEG(S));         % noise metrics on pop-free data
    assert(infoS.robustSd < 0.3 * noise, 'the virtual channel is much quieter than one channel (%.2f)', infoS.robustSd);
    assert(abs(infoS.medianChannelRobustSd - noise) < 0.2 * noise, 'per-channel robust sd ~ %g (got %.2f)', noise, infoS.medianChannelRobustSd);

    %% 3. per-channel DC offsets do not matter
    offsets = 1e5 * randn(N, 1);
    ref3 = RELAX_reference_channel(mkEEG(X + offsets));
    assert(max(abs(ref3.data - ref.data)) < 1e-6, 'offsets change nothing');

    %% 4. a minority of channels with their own pops does not move the median
    X4 = X; X4(1:10, 2000:end) = X4(1:10, 2000:end) + 2000;
    ref4 = RELAX_reference_channel(mkEEG(X4));
    assert(max(abs(demean(ref4.data) - demean(ref.data))) < 10);

    %% 5. demeaned per segment: offsets that differ between segments do not leak
    X5 = S; X5(:, 6001:end) = X5(:, 6001:end) + 1000 * (1:N)';
    [ref5, info5] = RELAX_reference_channel(mkEEG(X5, 6000.5));
    assert(isequal(info5.segments, [1 6000; 6001 n]));
    assert(max(abs(ref5.data)) < 10, 'no step at the boundary from the offsets (got %.1f)', max(abs(ref5.data)));

    %% 6. flat and explicitly excluded channels are left out
    X6 = X; X6(1, :) = 0; X6(2, 3000:end) = X6(2, 3000:end) + 1e4;
    [ref6, info6] = RELAX_reference_channel(mkEEG(X6), 'ExcludeLabels', {'E2'});
    assert(info6.nChannelsUsed == N - 2);
    assert(isequal(sort(info6.excludedLabels(:)), {'E1'; 'E2'}));
    assert(~any(ismember({'E1','E2'}, info6.usedLabels)));
    assert(max(abs(demean(ref6.data) - demean(ref.data))) < 10);

    %% 7. too few channels
    [ref7, info7] = RELAX_reference_channel(mkEEG(X(1:4, :)), 'MinChannels', 8);
    assert(~info7.available && ~isempty(info7.status));
    assert(all(isnan(ref7.data)) && isequal(size(ref7.data), [1 n]));

    %% 8. single precision input comes out double
    ref8 = RELAX_reference_channel(mkEEG(single(X)));
    assert(isa(ref8.data, 'double'));

    fprintf('test_RELAX_reference_channel: PASS\n');
end
