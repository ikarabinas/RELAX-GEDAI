% GEDAI_alpha_mask_psd() - Welch PSD from the runs of true in a sample mask.
%
% Usage:
%   >> [freqs, power] = GEDAI_alpha_mask_psd(data, srate, mask, winSec, maxHz);
%
% Each run of true in mask is one segment for RELAX_QC_segment_psd, so no window spans
% a gap and anything inside a gap (steps, excised data) cannot reach the spectrum.
%
% Outputs:
%   freqs - [1 x F]
%   power - [nChannels x F]

function [freqs, power] = GEDAI_alpha_mask_psd(data, srate, mask, winSec, maxHz)

    paddedMask = [false, mask(:)', false];
    segments = [find(diff(paddedMask) == 1)', find(diff(paddedMask) == -1)'];
    [freqs, power] = RELAX_QC_segment_psd(data, srate, segments, winSec, maxHz);
end
