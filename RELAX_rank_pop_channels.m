function info = RELAX_rank_pop_channels(steps, info, keep)
% RELAX_rank_pop_channels - Rank chronically bad contacts from a chosen subset of steps.
%
%   info = RELAX_rank_pop_channels(steps, info)
%   info = RELAX_rank_pop_channels(steps, info, keep)
%
% Recomputes info.stepRatePerMinute, info.transitionFraction, info.badChannels and
% info.badChannelReason (worst-first, as RELAX_channel_budget 'spend' requires) counting
% only the steps with keep(i) true as the channel's own. RELAX_detect_steps calls this with
% keep = ~steps.commonMode; the wrapper calls it again after RELAX_flag_eye_steps with the
% default keep = ~steps.commonMode & ~steps.eyeMovement, because a saccade is a genuine
% DC step on the periocular channels and would otherwise reject the very electrodes the
% blink / HEOG stages need (m159: E10 E46 E18 E37 ... all "bad" by rate).
%
% steps - RELAX_detect_steps output (.channel .onset .offset .commonMode [.eyeMovement])
% info  - RELAX_detect_steps info (.params .labels .minutes .pnts
%         .degenerateChannels - these are listed first, reason 'flat')
% keep  - logical N x 1, default ~commonMode & ~eyeMovement (missing fields = all false)

    n = numel(steps.channel);
    if nargin < 3 || isempty(keep)
        keep = true(n, 1);
        if isfield(steps, 'commonMode') && ~isempty(steps.commonMode)
            keep = keep & ~logical(steps.commonMode(:));
        end
        if isfield(steps, 'eyeMovement') && ~isempty(steps.eyeMovement)
            keep = keep & ~logical(steps.eyeMovement(:));
        end
    end
    keep = logical(keep(:));
    nb = numel(info.labels);

    ownCh  = reshape(steps.channel(keep), [], 1);
    ownLen = reshape(steps.offset(keep) - steps.onset(keep) + 1, [], 1);
    info.stepRatePerMinute  = accumarray(ownCh, 1, [nb 1]) / max(info.minutes, eps);
    info.transitionFraction = accumarray(ownCh, ownLen, [nb 1]) / info.pnts;

    R = info.params.StepsPerMinuteToRejectChannel;
    F = info.params.MaxTransitionFractionToRejectChannel;
    scores = [info.stepRatePerMinute / R, info.transitionFraction / F];
    [score, which] = max(scores, [], 2);
    % degenerate (flat / railed) channels first - unambiguous, and left in they crash
    % pop_jointprob downstream; then the scored channels worst-first
    flat = find(ismember(info.labels(:), info.degenerateChannels(:)));
    idx = find(score > 1 & ~ismember((1:nb)', flat));
    [~, ord] = sort(score(idx), 'descend');
    idx = [flat(:); idx(ord)];
    info.badChannels = reshape(info.labels(idx), [], 1);
    reason = cell(numel(idx), 1);
    for k = 1:numel(idx)
        c = idx(k);
        if k <= numel(flat)
            reason{k} = 'flat: constant within every segment';
            continue
        end
        switch which(c)
            case 1, reason{k} = sprintf('%.1f steps/min exceeds %g', info.stepRatePerMinute(c), R);
            case 2, reason{k} = sprintf('%.3f of samples in transition exceeds %g', info.transitionFraction(c), F);
        end
    end
    info.badChannelReason = reason;
end
