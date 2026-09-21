load('/home/nga4004/gitRepos/RELAX-GEDAI/reference_pop_survey.mat', 'results');
for k = 1:numel(results)
    r = results{k}; c = r.ref{1}.candidates; g = r.ref{1}.gated;
    fprintf('\n%s: %d candidates, %d gated, srate %g\n', r.file, numel(c.sample), sum(g), r.srate);
    key = [c.onset(:) c.offset(:)];
    [~, ~, j] = unique(key, 'rows'); dup = accumarray(j, 1) > 1;
    fprintf('  duplicate (onset,offset) pairs: %d\n', sum(dup));
    for d = find(dup)'
        idx = find(j == d)';
        fprintf('    cands %s: onset %d offset %d amp %s z %s gated %s\n', mat2str(idx), key(idx(1),1), key(idx(1),2), ...
            mat2str(round(c.amplitude(idx)',1)), mat2str(round(c.z(idx)',1)), mat2str(g(idx)'));
    end
    gi = find(g); 
    if numel(gi) > 1
        gaps = diff(c.sample(gi)) / r.srate;
        fprintf('  gated: %d within 2 s of the previous gated candidate (episodes: %d)\n', sum(gaps < 2), 1 + sum(gaps >= 2));
    end
    overl = 0;
    for a = 1:numel(gi), for b = a+1:numel(gi)
        if c.onset(gi(b)) <= c.offset(gi(a)) + 0 && c.offset(gi(b)) >= c.onset(gi(a)); overl = overl + 1; end
    end, end
    fprintf('  gated pairs with overlapping [onset,offset]: %d\n', overl);
    fl10 = r.ref{2}.candidates; g10 = r.ref{2}.gated;
    fprintf('  floor 10: %d gated, |amp| of gated: %s\n', sum(g10), mat2str(round(sort(abs(fl10.amplitude(g10)))',0)));
    u = r.unmatchedCommonMode;
    fprintf('  unmatched common-mode: %d; agree>=0.9&spread<=0.5: %d; median agree %.2f spread %.2f near %.2f\n', numel(u.lo), ...
        sum(u.signAgreement >= 0.9 & u.amplitudeSpread <= 0.5), median(u.signAgreement,'omitnan'), median(u.amplitudeSpread,'omitnan'), median(u.fractionNear,'omitnan'));
end
