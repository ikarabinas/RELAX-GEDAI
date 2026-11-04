% RELAX_IC_muscle_slopes() - Log-frequency log-power slope of each ICA component over
%                           12-70 Hz (excluding +/-6 Hz around line noise), from one
%                           Welch spectrum per component. Slopes at or above
%                           RELAX_cfg.MuscleSlopeThreshold are what targeted wICA
%                           treats as muscle.
%
% Usage:
%   >> slopes = RELAX_IC_muscle_slopes(icaact, srate, lineNoiseFrequency);
%   >> slopes = RELAX_IC_muscle_slopes(icaact, srate, lineNoiseFrequency, segments);
%
% Inputs:
%   icaact             - components x samples (x trials) activations
%   srate              - sampling rate in Hz
%   lineNoiseFrequency - RELAX_cfg.LineNoiseFrequency
%   segments           - (optional) [nSegments x 2] start/end samples of contiguous data, from
%                         RELAX_get_segments. Each segment of at least two seconds gets its
%                         own spectrum and the spectra are averaged weighted by length, so a
%                         step where two segments meet adds no broadband power to the fit.
%                         Without it the record is treated as one contiguous signal.
%
% Outputs:
%   slopes  - 1 x components
%   fftBins - (optional) components x numel(freq) power in 0.5 Hz bins that the fit used
%   freq    - (optional) 1 x numel(freq) bin centres in Hz

function [slopes, fftBins, freq] = RELAX_IC_muscle_slopes(icaact, srate, lineNoiseFrequency, segments)

    minimumSegmentSeconds = 2;
    muscleFreqIn = [12,70];
    freq = 1:0.5:100;

    eegData = reshape(icaact,size(icaact,1),[]);
    if nargin < 4
        segments = [1, size(eegData,2)];
    end
    segmentLengths = segments(:,2) - segments(:,1) + 1;
    usable = segmentLengths >= minimumSegmentSeconds * srate;
    if ~any(usable)
        error('RELAX:noSegmentLongEnough', 'No segment is at least %g s long, so no spectrum can be estimated.', minimumSegmentSeconds);
    end
    segments = segments(usable,:);
    segmentLengths = segmentLengths(usable);

    fftBins = zeros(size(eegData,1), numel(freq));
    for s = 1:size(segments,1)
        segmentBins = binnedSpectrum(eegData(:,segments(s,1):segments(s,2)), srate, freq);
        fftBins = fftBins + segmentBins * segmentLengths(s) / sum(segmentLengths);
    end

    muscleFreqEx = [lineNoiseFrequency-6 lineNoiseFrequency+6];
    [~,fin1] = min(abs(muscleFreqIn(1) - freq));
    [~,fin2] = min(abs(muscleFreqIn(2) - freq));
    freqHz = freq(fin1:fin2);
    [~,fex1] = min(abs(muscleFreqEx(1) - freqHz));
    [~,fex2] = min(abs(muscleFreqEx(2) - freqHz));
    freqHz(fex1:fex2) = [];

    comps = size(eegData,1);
    slopes = zeros(1, comps);
    for compNum = 1:comps
        freqPow = fftBins(compNum,fin1:fin2);
        freqPow(fex1:fex2) = [];
        p = polyfit(log(freqHz),log(freqPow),1);
        slopes(compNum) = p(1);
    end
end

% Power in 0.5 Hz bins centred on each frequency in freq, from one Welch spectrum whose window
% is the whole of the data given
function fftBins = binnedSpectrum(eegData, srate, freq)
    [pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),srate);
    FFTout = pxx';
    fp = fp';
    fftBins = zeros(size(FFTout,1),numel(freq));
    for a = 1:numel(freq)
        [~, index1] = min(abs(fp-(freq(a)-0.25)));
        [~, index2] = min(abs(fp-(freq(a)+0.25)));
        fftBins(:,a) = mean(FFTout(:,index1:index2),2);
    end
end
