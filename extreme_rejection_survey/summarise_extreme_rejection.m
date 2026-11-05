function summarise_extreme_rejection(resultsDir)
% Folds the per-file survey results into <resultsDir>/../extreme_rejection_survey.csv and
% prints the across-file and within-participant comparison of rejection rules.

    listing = dir(fullfile(resultsDir, '*.mat'));
    rows = cell(numel(listing), 1);
    badChannelCounts = [];
    for i = 1:numel(listing)
        loaded = load(fullfile(listing(i).folder, listing(i).name), 'result');
        rows{i} = file_row(loaded.result);
        flaggedEpochs = loaded.result.badChannelsPerEpoch(loaded.result.badChannelsPerEpoch > 0);
        badChannelCounts = [badChannelCounts; double(flaggedEpochs(:))]; %#ok<AGROW>
    end
    T = struct2table([rows{:}]');
    T = sortrows(T, 'file');
    outCsv = fullfile(fileparts(resultsDir), 'extreme_rejection_survey.csv');
    writetable(T, outCsv);
    fprintf('Wrote %s (%d files)\n\n', outCsv, height(T));

    fprintf('== Reproduction check ==\n');
    fprintf('Files whose rerun matches the saved proportion excluded (|diff| < 1e-6): %d / %d\n', ...
        sum(abs(T.rerunExcluded - T.savedExcluded) < 1e-6), height(T));
    fprintf('Largest |rerun - saved| proportion excluded: %.2g\n', max(abs(T.rerunExcluded - T.savedExcluded)));
    fprintf('Files where per-channel flags disagree with RELAX epoch marks: %d\n', sum(T.epochMismatches > 0));
    fprintf('Files where the run''s own rule differs from the rerun: %d\n\n', sum(~T.ruleMatchesRerun));

    fprintf('== Proportion of evaluated data excluded, by rule (median [IQR] across files; mean within-participant SD across sessions) ==\n');
    ruleCols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'excluded_'));
    for c = ruleCols
        x = T.(c{1});
        fprintf('%-22s %.3f [%.3f-%.3f]  max %.3f  within-ppt SD %.3f\n', erase(c{1}, 'excluded_'), ...
            median(x), prctile(x, 25), prctile(x, 75), max(x), within_participant_sd(T.participant, x));
    end
    fprintf('Total deleted incl. never-evaluated data, current rule: median %.3f, max %.3f\n\n', ...
        median(T.deleted_k1), max(T.deleted_k1));

    fprintf('== Per criterion (median across files) ==\n');
    critCols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'flagged_'));
    for c = critCols
        name = erase(c{1}, 'flagged_');
        fprintf('%-12s epochs flagged %.3f   flagged by this criterion alone %.3f\n', name, ...
            median(T.(c{1})), median(T.(['only_' name])));
    end
    fprintf('\n== Bad channels per flagged epoch, pooled over files ==\n');
    edges = [1 2 3 5 10 20 inf];
    for e = 1:numel(edges) - 1
        fprintf('%3d-%-4s channels: %.3f of flagged epochs\n', edges(e), num2str(edges(e + 1) - 1), ...
            mean(badChannelCounts >= edges(e) & badChannelCounts < edges(e + 1)));
    end
    fprintf('\nCorrelation of current-rule loss with channel count: r = %.2f\n', corr(T.nChan, T.excluded_k1));

    fprintf('\n== TMS analysis windows (post 0.5-2.5 s, pre -2.1 to -0.1 s) ==\n');
    fprintf('Files using trainIdx boundaries: %d / %d (others fall back to all boundaries)\n', sum(T.usedTrainIdx), height(T));
    fprintf('Windows: %d post, %d pre (median per file %g / %g); dropped for crossing another boundary: %d post, %d pre\n', ...
        sum(T.nPost), sum(T.nPre), median(T.nPost), median(T.nPre), sum(T.nDroppedPost), sum(T.nDroppedPre));
    fprintf('Touching never-evaluated data (lost under every rule, not counted below): %d post, %d pre\n\n', ...
        sum(T.nPostNeverEvaluated), sum(T.nPreNeverEvaluated));

    fprintf('== Rejected stretches and windows losing data, by rule (never-evaluated data excluded) ==\n');
    fprintf('%-12s %10s %10s %10s %9s %16s %16s %14s\n', 'rule', 'stretches', 'per file', 'splitting', 'per file', ...
        'post hit (%)', 'pre hit (%)', 'worst file hit');
    for c = ruleCols
        name = erase(c{1}, 'excluded_');
        post = T.(['postHit_' name]);
        pre = T.(['preHit_' name]);
        worst = max((post + pre) ./ (T.nPost + T.nPre));
        fprintf('%-12s %10d %10g %10d %9g %8d (%5.2f) %8d (%5.2f) %13.1f%%\n', name, ...
            sum(T.(['stretches_' name])), median(T.(['stretches_' name])), ...
            sum(T.(['splitting_' name])), median(T.(['splitting_' name])), ...
            sum(post), 100 * sum(post) / sum(T.nPost), sum(pre), 100 * sum(pre) / sum(T.nPre), 100 * worst);
    end
    fprintf('(per file = median across files; worst file hit = largest share of a file''s post+pre windows losing data)\n');
end

function row = file_row(r)
    row.file = r.file;
    row.participant = extractBefore(r.file, '_');
    row.nChan = r.nChan;
    row.minutes = r.minutes;
    row.savedExcluded = r.check.savedProportionExcluded;
    row.rerunExcluded = r.check.rerunProportionExcluded;
    row.savedDeleted = r.check.savedProportionDeleted;
    row.epochMismatches = sum(cell2mat(struct2cell(r.check.mismatchedEpochs)));
    row.ruleMatchesRerun = r.check.ruleMatchesRerun;
    for i = 1:numel(r.rule)
        name = rule_field(r.rule(i).name);
        row.(['excluded_' name]) = r.rule(i).proportionExcluded;
        row = add_geometry(row, name, r.rule(i).geometry);
    end
    w = r.windows;
    row.usedTrainIdx = w.usedTrainIdx;
    row.nPost = w.nPost;
    row.nPre = w.nPre;
    row.nDroppedPost = w.nDroppedPost;
    row.nDroppedPre = w.nDroppedPre;
    row.nPostNeverEvaluated = w.nPostNeverEvaluated;
    row.nPreNeverEvaluated = w.nPreNeverEvaluated;
    row.deleted_k1 = r.rule(1).proportionDeleted;
    for c = fieldnames(r.byCriterion)'
        row.(['flagged_' c{1}]) = r.byCriterion.(c{1}).proportionEpochsFlagged;
        row.(['only_' c{1}]) = r.byCriterion.(c{1}).proportionEpochsFlaggedOnlyByThis;
    end
    for d = r.dropWorst
        name = sprintf('dropWorst%d', d.nChannelsDropped);
        row.(['excluded_' name]) = d.proportionExcluded;
        row = add_geometry(row, name, d.geometry);
    end
    row.medianFlagsPerChannel = median(r.flagsPerChannel);
    row.maxFlagsPerChannel = max(r.flagsPerChannel);
end

function row = add_geometry(row, name, g)
    row.(['stretches_' name]) = g.nStretches;
    row.(['splitting_' name]) = g.nSplitting;
    row.(['postHit_' name]) = g.nPostHit;
    row.(['preHit_' name]) = g.nPreHit;
end

function name = rule_field(ruleName)
% 'k>=5' -> 'k5'; '>=2%(k=4)' -> 'pct2'. The pct rules' k varies per file, so they are
% named by percentage.
    if startsWith(ruleName, 'k')
        name = ['k' extractAfter(ruleName, '>=')];
    else
        name = ['pct' extractBetween(ruleName, '>=', '%')];
        name = [name{:}];
    end
end

function s = within_participant_sd(participant, x)
    [groups, ~] = findgroups(cellstr(participant));  % a one-row table holds it as char
    sds = splitapply(@(v) std(v), x, groups);
    counts = splitapply(@numel, x, groups);
    s = mean(sds(counts > 1));
end
