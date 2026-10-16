% RELAX_QC_pick_examples() - Choose the files shown as Figure 1 examples.
%
% Usage:
%   >> picks = RELAX_QC_pick_examples(files, scores, overrides);
%
% Automatic choice: the median file (its mean rank across the three scores is closest
% to the median of all files' mean ranks) and the worst file by blink amplitude, muscle
% power and pop count. A file filling several roles appears once with the roles joined.
% A non-empty overrides list replaces the automatic choice entirely; an override naming
% a file not in `files` warns (RELAX_QC:unknownOverride) and is dropped.

function picks = RELAX_QC_pick_examples(files, scores, overrides)

    picks = struct('file', {}, 'role', {}, 'centres', {});

    if ~isempty(overrides)
        for o = overrides(:)'
            file = string(o.file);
            if ~any(files == file)
                warning('RELAX_QC:unknownOverride', 'Override names %s, which is not in the dataset; ignored.', file);
                continue
            end
            k = find([string.empty, picks.file] == file, 1);
            if isempty(k)
                picks(end+1) = newPick(file, "chosen"); %#ok<AGROW>
                k = numel(picks);
            end
            picks(k).centres.(char(o.kind)) = o.centreSeconds;
        end
        return
    end

    meanRank = mean([rankOf(scores.blinkAmplitude), rankOf(scores.musclePower), rankOf(scores.popSteps)], 2, 'omitnan');
    if all(isnan(meanRank)), return, end
    [~, iMedian] = min(abs(meanRank - median(meanRank, 'omitnan')));
    candidates = {iMedian, "median"; worst(scores.blinkAmplitude), "worst blink"; ...
                  worst(scores.musclePower), "worst muscle"; worst(scores.popSteps), "most pops"};
    for c = 1:size(candidates, 1)
        i = candidates{c, 1};
        if isempty(i), continue, end
        k = find([string.empty, picks.file] == files(i), 1);
        if isempty(k)
            picks(end+1) = newPick(files(i), candidates{c, 2}); %#ok<AGROW>
        else
            picks(k).role = picks(k).role + ", " + candidates{c, 2};
        end
    end
end

function p = newPick(file, role)
    p = struct('file', file, 'role', role, 'centres', struct('blink', NaN, 'muscle', NaN, 'pop', NaN));
end

function r = rankOf(x)
    r = nan(size(x));
    scored = find(~isnan(x));
    [~, order] = sort(x(scored));
    r(scored(order)) = 1:numel(scored);
end

function i = worst(x)
    i = [];
    if any(~isnan(x)), [~, i] = max(x); end
end
