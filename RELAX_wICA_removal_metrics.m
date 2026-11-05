% RELAX_wICA_removal_metrics() - How much of each frequency band the wICA step removed, split
%                                into clean, muscle and blink data.
%
% Usage:
%   >> [metrics, perChannel] = RELAX_wICA_removal_metrics(EEGbefore, EEGafter, RELAX_cfg, blinkMask);
%
% Inputs:
%   EEGbefore - continuous data going into wICA
%   EEGafter  - the same data coming out of it (same channels and samples)
%   RELAX_cfg - uses MuscleSlopeThreshold (-0.31 when absent)
%   blinkMask - [1 x pnts] logical, 1 inside blinks (the wrapper passes RELAX_blinks_IQR_method's
%               eyeblinkmask on the wICA input); all zeros when omitted
%
% The removed signal is EEGbefore - EEGafter. Each stretch between boundary events is cut into
% non-overlapping 1 s epochs (so no epoch spans a boundary) and Hann-windowed periodograms are
% taken of the input and of the removed signal. Every epoch touching a blink is Blink. In the
% others each channel is judged on its own: Muscle where that channel's 7-75 Hz log-log input
% slope reaches MuscleSlopeThreshold, Clean otherwise. Judging per channel matters on a dense net:
% with 200 channels almost every epoch has some channel over the threshold. Brain signal lost to
% wICA shows up as removed power in the Clean data, chiefly in the alpha and beta bands.
%
% Outputs:
%   metrics    - scalar struct: RemovedFraction_<Clean|Muscle|Blink>_<band> = removed power /
%                input power summed over the channel-epochs of that kind (NaN when there are
%                none), plus NumberOfEpochs, ProportionEpochsBlink and
%                ProportionChannelEpochsMuscle (among epochs without blinks)
%   perChannel - labels, chanlocs, bandNames, bandEdgesHz, and [channels x bands] removed and
%                input power for each kind, for topographies

function [metrics, perChannel] = RELAX_wICA_removal_metrics(EEGbefore, EEGafter, RELAX_cfg, blinkMask)

    muscleSlopeThreshold = -0.31;
    if isfield(RELAX_cfg, 'MuscleSlopeThreshold'), muscleSlopeThreshold = RELAX_cfg.MuscleSlopeThreshold; end
    if nargin < 4 || isempty(blinkMask), blinkMask = false(1, EEGbefore.pnts); end

    bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'HighGamma', 'Broadband'};
    bandEdgesHz = [1 4; 4 8; 8 13; 13 30; 30 55; 65 100; 1 100];
    kindNames = {'Clean', 'Muscle', 'Blink'};

    epochLength = round(EEGbefore.srate);
    freqHz = (0:floor(epochLength/2)) * EEGbefore.srate / epochLength;
    taper = hann(epochLength)';
    inBand = (freqHz >= bandEdgesHz(:,1)) & (freqHz < bandEdgesHz(:,2));   % [bands x bins]
    slopeBins = freqHz >= 7 & freqHz <= 75;
    logSlopeFreq = log(freqHz(slopeBins));

    nChannels = EEGbefore.nbchan;
    removedPower = zeros(nChannels, numel(bandNames), numel(kindNames));
    inputPower = zeros(nChannels, numel(bandNames), numel(kindNames));
    channelEpochs = zeros(nChannels, numel(kindNames));
    nEpochs = 0;
    nBlinkEpochs = 0;

    segments = RELAX_get_segments(EEGbefore);
    for s = 1:size(segments, 1)
        for first = segments(s,1) : epochLength : segments(s,2) - epochLength + 1
            samples = first : first + epochLength - 1;
            input = double(EEGbefore.data(:, samples));
            removed = input - double(EEGafter.data(:, samples));
            inputBands = periodogram(input, taper, numel(freqHz));
            removedBands = periodogram(removed, taper, numel(freqHz)) * inBand';
            nEpochs = nEpochs + 1;

            if any(blinkMask(samples))
                kindOfChannel = 3 * ones(nChannels, 1);
                nBlinkEpochs = nBlinkEpochs + 1;
            else
                slopes = loglogSlopes(logSlopeFreq, log(inputBands(:, slopeBins)));
                kindOfChannel = 1 + (slopes >= muscleSlopeThreshold);
            end
            inputBands = inputBands * inBand';
            for kind = 1:numel(kindNames)
                rows = kindOfChannel == kind;
                removedPower(rows,:,kind) = removedPower(rows,:,kind) + removedBands(rows,:);
                inputPower(rows,:,kind) = inputPower(rows,:,kind) + inputBands(rows,:);
                channelEpochs(:,kind) = channelEpochs(:,kind) + rows;
            end
        end
    end

    for kind = 1:numel(kindNames)
        fraction = sum(removedPower(:,:,kind), 1) ./ sum(inputPower(:,:,kind), 1);
        if ~any(channelEpochs(:,kind)), fraction(:) = NaN; end
        for b = 1:numel(bandNames)
            metrics.(['RemovedFraction_' kindNames{kind} '_' bandNames{b}]) = fraction(b);
        end
    end
    metrics.NumberOfEpochs = nEpochs;
    metrics.ProportionEpochsBlink = nBlinkEpochs / max(nEpochs, 1);
    metrics.ProportionChannelEpochsMuscle = sum(channelEpochs(:,2)) / max(sum(channelEpochs(:,1:2), 'all'), 1);

    perChannel.labels = {EEGbefore.chanlocs.labels};
    perChannel.chanlocs = EEGbefore.chanlocs;
    perChannel.bandNames = bandNames;
    perChannel.bandEdgesHz = bandEdgesHz;
    for kind = 1:numel(kindNames)
        perChannel.(['channelEpochs' kindNames{kind}]) = channelEpochs(:,kind);
        perChannel.(['removedPower' kindNames{kind}]) = removedPower(:,:,kind);
        perChannel.(['inputPower' kindNames{kind}]) = inputPower(:,:,kind);
    end
end

function spectrum = periodogram(data, taper, nBins)
    data = data - mean(data, 2);
    spectrum = abs(fft(data .* taper, [], 2)).^2;
    spectrum = spectrum(:, 1:nBins);
end

function slopes = loglogSlopes(logFreq, logPower)
    centredFreq = logFreq - mean(logFreq);
    slopes = (logPower - mean(logPower, 2)) * centredFreq' / sum(centredFreq.^2);
end
