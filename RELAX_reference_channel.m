% RELAX_reference_channel() - Build the virtual reference channel: the cross-channel median
%                             of per-segment-demeaned data.
%
% Usage:
%   >> [ref, info] = RELAX_reference_channel(EEG);
%   >> [ref, info] = RELAX_reference_channel(EEG, 'ExcludeLabels', {'E2'}, 'MinChannels', 8);
%
% On a Cz-referenced montage every channel is x_i(t) = s_i(t) - r(t), so a pop on the
% reference contact is the same waveform on every channel. The median across channels
% keeps that shared waveform and suppresses everything channel-specific: with ~200
% channels its noise is ~1/sqrt(N) of one channel's, and a minority of channels with their
% own concurrent artifact does not move it. Demeaning is per inter-boundary segment
% because raw DC-coupled channels sit at offsets tens of mV apart; without it the median
% would hop between channels as they reorder.
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   Name-value pairs (optional):
%     ExcludeLabels - channel labels left out of the median (default {})
%     MinChannels   - fewer usable channels than this -> info.available = false (default 8)
%
% Outputs:
%   ref  - 1-channel EEG struct (.data 1 x pnts double, .nbchan 1, .pnts, .srate, .trials,
%          .xmin, .xmax, .event copied from EEG so boundaries are preserved, .chanlocs with
%          label 'REFMEDIAN'); NaN data when unavailable
%   info - .available, .status, .nChannelsUsed, .usedLabels, .excludedLabels (caller's
%          plus flat channels, RELAX_flat_channels), .segments, .robustSd (of ref.data,
%          segment-length weighted), .medianChannelRobustSd (median over used channels of
%          their own demeaned robust sd - the noise the median is averaging over)

function [ref, info] = RELAX_reference_channel(EEG, varargin)
    p = inputParser;
    p.addParameter('ExcludeLabels', {});
    p.addParameter('MinChannels', 8);
    p.parse(varargin{:});
    opt = p.Results;

    n = EEG.pnts;
    labels = reshape({EEG.chanlocs.labels}, [], 1);
    segments = RELAX_get_segments(EEG);
    excluded = ismember(labels, cellstr(opt.ExcludeLabels));
    flat = RELAX_flat_channels(EEG.data, segments);
    used = find(~excluded & ~flat);

    ref = struct('data', nan(1, n), 'nbchan', 1, 'pnts', n, 'srate', EEG.srate, ...
                 'trials', 1, 'xmin', 0, 'xmax', (n - 1) / EEG.srate, ...
                 'event', EEG.event, 'chanlocs', struct('labels', 'REFMEDIAN'));
    if isfield(EEG, 'xmin'); ref.xmin = EEG.xmin; end
    if isfield(EEG, 'xmax'); ref.xmax = EEG.xmax; end
    info = struct('available', false, 'status', '', 'nChannelsUsed', numel(used), ...
                  'usedLabels', {labels(used)}, 'excludedLabels', {labels(excluded | flat)}, ...
                  'segments', segments, 'robustSd', NaN, 'medianChannelRobustSd', NaN);
    if numel(used) < opt.MinChannels
        info.status = sprintf('%d usable channel(s), fewer than MinChannels = %d', ...
                              numel(used), opt.MinChannels);
        return
    end

    nSeg = size(segments, 1);
    segLen = (segments(:,2) - segments(:,1) + 1)';
    channelSd = nan(numel(used), nSeg);
    refSd = nan(1, nSeg);
    for s = 1:nSeg
        a = segments(s,1); b = segments(s,2);
        X = double(EEG.data(used, a:b));
        X = X - median(X, 2, 'omitnan');
        ref.data(a:b) = median(X, 1, 'omitnan');
        channelSd(:, s) = robustSdRows(X);
        refSd(s) = robustSdRows(ref.data(a:b));
    end
    info.robustSd = sum(refSd .* segLen) / sum(segLen);
    info.medianChannelRobustSd = median(channelSd * segLen' / sum(segLen));
    info.available = true;
    info.status = 'ok';
end

function s = robustSdRows(X)
    s = 1.4826 * median(abs(X - median(X, 2, 'omitnan')), 2, 'omitnan');
end
