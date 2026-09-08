% RELAX_pad_blink_mask_block() - Dilate and smooth a binary blink mask over one block.
%
% Usage:
%   >> out = RELAX_pad_blink_mask_block(mask, moving_mean_length);
%
% Inputs:
%   mask               - [1 x nSamples] binary mask of blink periods
%   moving_mean_length - dilation radius in samples
%
% Outputs:
%   out - [1 x nSamples] mask dilated by moving_mean_length in each direction and then
%         smoothed with a centred moving average of the same half-width, giving the
%         tapered weighting that wICA multiplies blink components by.
%
% Applied one block at a time so that, when called through RELAX_apply_by_segments, a
% blink near the end of one segment cannot dilate into the next.

function out = RELAX_pad_blink_mask_block(mask, moving_mean_length)

    padded_blink_periods = double(mask);
    n = size(padded_blink_periods,2);

    for c = flip(1:max(0, n-(moving_mean_length+1)))
        if padded_blink_periods(1,c)==1
            padded_blink_periods(1,c:c+moving_mean_length)=1;
        end
    end
    for c = (moving_mean_length+1):n
        if padded_blink_periods(1,c)==1
            padded_blink_periods(1,c-moving_mean_length:c)=1;
        end
    end

    out = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
end
