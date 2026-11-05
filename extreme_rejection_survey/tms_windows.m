function windows = tms_windows(EEG, segments)
% Sample ranges of the downstream analysis epochs: 0.5 to 2.5 s after each TMS
% boundary (post) and 2.1 to 0.1 s before it (pre). TMS boundaries are the 'boundary'
% events carrying a trainIdx (the ones RELAX_Wrapper_slurm.m turns into TMS_train_N
% markers); if none carry one, every boundary is used and usedTrainIdx says so.
% A window reaching across another boundary is not a real analysis epoch, so it is
% dropped and counted in nDroppedPost/nDroppedPre.

    isBoundary = arrayfun(@(x) (ischar(x.type) || isstring(x.type)) && strcmp(x.type, 'boundary'), EEG.event);
    isTms = isBoundary & arrayfun(@(x) isfield(x, 'trainIdx') && ~isempty(x.trainIdx), EEG.event);
    windows.usedTrainIdx = any(isTms);
    if ~windows.usedTrainIdx
        isTms = isBoundary;
    end
    b = sort([EEG.event(isTms).latency]);
    b = b(b > 1 & b < EEG.pnts);
    firstAfter = floor(b)' + 1;   % same convention as RELAX_get_segments
    lastBefore = floor(b)';

    post = [firstAfter + round(0.5 * EEG.srate), firstAfter + round(2.5 * EEG.srate) - 1];
    pre  = [lastBefore - round(2.1 * EEG.srate) + 1, lastBefore - round(0.1 * EEG.srate)];
    keepPost = in_one_segment(post, firstAfter, segments);
    keepPre  = in_one_segment(pre, lastBefore, segments);
    windows.post = post(keepPost, :);
    windows.pre  = pre(keepPre, :);
    windows.nDroppedPost = sum(~keepPost);
    windows.nDroppedPre  = sum(~keepPre);
end

function ok = in_one_segment(ranges, anchor, segments)
    ok = false(size(ranges, 1), 1);
    for i = 1:size(ranges, 1)
        s = find(segments(:,1) <= anchor(i) & segments(:,2) >= anchor(i), 1);
        ok(i) = ~isempty(s) && ranges(i,1) >= segments(s,1) && ranges(i,2) <= segments(s,2);
    end
end
