% RELAX_step_uniformity() - Level change across one interval on every channel, and how
%                           uniform it is across channels.
%
% Usage:
%   >> [amp, summary] = RELAX_step_uniformity(EEG.data, onset, offset, W, segStart, segEnd);
%
% A reference-electrode pop is the same level change on every channel (x_i = s_i - r); a
% head-movement event is not. amp(i) is median(x_i(offset : offset+W-1)) minus
% median(x_i(onset-W+1 : onset)), both windows clipped to the segment.
%
% Inputs:
%   data           - nbchan x pnts
%   onset, offset  - last pre-level sample and first post-level sample (RELAX_detect_steps)
%   W              - level window in samples
%   segStart, segEnd - the inter-boundary segment holding the interval (inclusive)
%
% Outputs:
%   amp     - nbchan x 1 double, NaN for a channel with no finite samples in a window
%   summary - .medianAmplitude, .nChannels (finite amps), .signAgreement (fraction whose
%             sign matches the median's), .amplitudeSpread (1.4826*MAD / |median|; Inf
%             when the median is 0), .fractionNear (fraction within 25 % of the median).
%             All NaN when no channel is finite.

function [amp, summary] = RELAX_step_uniformity(data, onset, offset, W, segStart, segEnd)
    preWin  = max(segStart, onset - W + 1) : onset;
    postWin = offset : min(segEnd, offset + W - 1);
    pre  = median(double(data(:, preWin)),  2, 'omitnan');
    post = median(double(data(:, postWin)), 2, 'omitnan');
    amp = post - pre;

    ok = isfinite(amp);
    summary = struct('medianAmplitude', NaN, 'nChannels', sum(ok), 'signAgreement', NaN, ...
                     'amplitudeSpread', NaN, 'fractionNear', NaN);
    if ~any(ok); return; end
    a = amp(ok);
    A = median(a);
    summary.medianAmplitude = A;
    summary.signAgreement   = mean(sign(a) == sign(A));
    summary.amplitudeSpread = 1.4826 * median(abs(a - A)) / abs(A);
    summary.fractionNear    = mean(abs(a - A) <= 0.25 * abs(A));
end
