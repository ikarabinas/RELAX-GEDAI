% RELAX_wavelet_denoise_segment() - Level-5 coif5 stationary wavelet thresholding of
%                                   one contiguous block of data.
%
% Usage:
%   >> comp_out = RELAX_wavelet_denoise_segment(comp_in, wavelet_threshold, threshold_type);
%
% Inputs:
%   comp_in           - [1 x nSamples] row vector (one independent component)
%   wavelet_threshold - scalar threshold, as returned by ddencmp
%   threshold_type    - 's' (soft) or 'h' (hard), as returned by ddencmp
%
% Outputs:
%   comp_out - [1 x nSamples], the same length as comp_in
%
% swt requires the signal length to be divisible by 2^5, so the block is zero-padded
% up to the next multiple of 32 and trimmed back afterwards. The threshold is passed
% in rather than computed here so that a caller processing several segments of the
% same component can apply one threshold derived from the whole component, instead of
% a different threshold per segment.

function comp_out = RELAX_wavelet_denoise_segment(comp_in, wavelet_threshold, threshold_type)

    check_padding_required = mod(numel(comp_in), 2^5);
    if check_padding_required ~= 0
        padding = zeros(1,(2^5)-check_padding_required);
    else
        padding = [];
    end

    padded_comp = [comp_in, padding];

    wavelet_transform = swt(padded_comp,5,'coif5');
    thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold);
    comp_out = iswt(thresholded_wavelet_transform,'coif5');

    if ~isempty(padding)
        comp_out = comp_out(1:end-numel(padding));
    end
end
