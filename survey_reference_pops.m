function results = survey_reference_pops(outFile, opts)
% SURVEY_REFERENCE_POPS - How common are reference-electrode pops on real files, and do
% they look like the model says (same sign, same amplitude, same waveform on every channel)?
%
% For each file, on the production 204-channel montage:
%   1. RELAX_detect_reference_pops with the uniformity gate OFF, at two amplitude floors
%      (opts.floors, uV), so every candidate on the cross-channel median is kept with its
%      signAgreement / amplitudeSpread / fractionNear diagnostics; the summary then reports
%      how many pass the shipped gate (0.9 / 0.5).
%   2. The existing per-channel detector (RELAX_detect_steps at the shipped config) so its
%      common-mode events can be matched against the reference candidates: what fraction
%      of common-mode events are reference pops, and what fraction of reference pops the
%      per-channel pass also saw.
%   3. Uniformity of every UNMATCHED common-mode event, to see whether a same-sign /
%      variable-amplitude class (ground pop) exists.
%   4. Figures: per file, up to opts.maxFigures candidates (largest |amplitude| first),
%      each showing the virtual channel and 12 random channels over [onset-1 s, offset+1 s]
%      and the sorted per-channel amplitudes; plus one scatter of amplitudeSpread vs
%      signAgreement over all files for reference candidates and unmatched common-mode
%      events.
%   5. If RELAX_subtract_reference_pops exists on the path (Task 5 of the plan), the
%      subtraction arm: apply it with the gated candidates at the first floor, re-run the
%      reference detector and the per-channel detector on the corrected data, and report
%      residual reference candidates and common-mode events before/after.
%
% Run on a compute node only: sbatch run_survey_reference_pops.sh
%
% opts (optional): .dataDir, .files, .montageFile, .floors ([30 10]), .maxFigures (6),
%   .figDir, .commonModeFraction (0.15), .slackMs (50), .subtract (true).

    if nargin < 1 || isempty(outFile); outFile = 'reference_pop_survey.mat'; end
    if nargin < 2; opts = struct(); end
    opts = setDefault(opts, 'dataDir', '~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Sep26/rest_treat_joined');
    opts = setDefault(opts, 'files', {'m159_dlpfc_day1_stitched.set', 'm497_dlpfc_day1_stitched.set', ...
        'm208_dlpfc_day1_stitched.set', 'm583_dlpfc_day1_stitched.set', 'm155_dlpfc_day1_stitched.set'});
    opts = setDefault(opts, 'montageFile', '/home/nga4004/TMS_response_variability/montage_204_egi_indices.txt');
    opts = setDefault(opts, 'floors', [30 10]);
    opts = setDefault(opts, 'maxFigures', 6);
    opts = setDefault(opts, 'figDir', fullfile('docs', 'superpowers', 'reports', 'reference_pop_survey'));
    opts = setDefault(opts, 'commonModeFraction', 0.15);
    opts = setDefault(opts, 'slackMs', 50);
    opts = setDefault(opts, 'subtract', true);
    if ~exist(opts.figDir, 'dir'); mkdir(opts.figDir); end

    eeglab nogui;
    RELAX_add_vendored_noisetools();
    electrodesToDelete = montageDeletions(opts.montageFile);
    detectOpts = {'WindowMs', 150, 'TransitionMs', 100, 'ZThreshold', 4, 'ChannelZThreshold', 4, ...
                  'RefractoryMs', 20, 'MaxTransitionMs', 300, 'MaxPlateauDriftFraction', Inf};
    canSubtract = opts.subtract && exist('RELAX_subtract_reference_pops', 'file') == 2;

    results = {};
    for k = 1:numel(opts.files)
        fprintf('\n===== %s\n', opts.files{k});
        EEG = pop_loadset('filename', opts.files{k}, 'filepath', opts.dataDir);
        EEG = pop_select(EEG, 'nochannel', intersect(electrodesToDelete, {EEG.chanlocs.labels}));
        EEG = eeg_checkset(EEG);
        r = struct('file', opts.files{k}, 'nbchan', EEG.nbchan, 'srate', EEG.srate, ...
                   'minutes', EEG.pnts / EEG.srate / 60, 'nSegments', size(RELAX_get_segments(EEG), 1));
        slack = round(opts.slackMs * EEG.srate / 1000);

        % 1. reference candidates, gate off, per floor
        r.floors = opts.floors;
        r.ref = cell(1, numel(opts.floors));
        for f = 1:numel(opts.floors)
            t0 = tic;
            [~, info] = RELAX_detect_reference_pops(EEG, 'MinSignAgreement', 0, 'MaxAmplitudeSpread', Inf, ...
                'Detect', [detectOpts, {'MinAmplitudeMicrovolts', opts.floors(f)}]);
            info.seconds = toc(t0);
            r.ref{f} = rmfield(info, 'refChannel');
            r.ref{f}.refTrace = single(info.refChannel.data);
            r.ref{f}.gated = info.candidates.signAgreement >= 0.9 & info.candidates.amplitudeSpread <= 0.5;
            % gate off -> every candidate accepted -> this holds all of them, for the figures
            r.ref{f}.perChannelAmplitudeAll = single(info.perChannelAmplitude);
            fprintf(['Reference channel (floor %g uV): %d usable channels, robust sd %.2f uV (channels %.2f); ' ...
                     '%d candidate(s) in %.1f min, %d pass the 0.9/0.5 gate; %.1f s\n'], opts.floors(f), ...
                info.channelInfo.nChannelsUsed, info.channelInfo.robustSd, info.channelInfo.medianChannelRobustSd, ...
                info.nCandidates, r.minutes, sum(r.ref{f}.gated), info.seconds);
            printCandidates(info.candidates, EEG.srate);
        end

        % 2. the per-channel detector at the shipped config, and the match
        t0 = tic;
        [steps, pinfo] = RELAX_detect_steps(EEG, detectOpts{:}, 'MinAmplitudeMicrovolts', 30, ...
            'CommonModeChannelFraction', opts.commonModeFraction);
        r.perChannel = struct('seconds', toc(t0), 'nSteps', pinfo.nSteps, 'events', pinfo.events, ...
                              'nCommonModeEvents', pinfo.nCommonModeEvents);
        cm = pinfo.events.commonMode;
        c1 = r.ref{1}.candidates;
        [refHasCM, cmHasRef] = RELAX_match_intervals(c1.onset, c1.offset, pinfo.events.lo(cm), pinfo.events.hi(cm), slack);
        r.match = struct('refHasCommonMode', refHasCM, 'commonModeHasRef', cmHasRef);
        fprintf('Per-channel detector: %d steps, %d events, %d common-mode; %d/%d reference candidates matched a common-mode event; %d/%d common-mode events matched a reference candidate (%.1f s)\n', ...
            pinfo.nSteps, numel(pinfo.events.id), pinfo.nCommonModeEvents, sum(refHasCM), numel(refHasCM), ...
            sum(cmHasRef), numel(cmHasRef), r.perChannel.seconds);

        % 3. uniformity of unmatched common-mode events
        W = pinfo.params.derived.W; segments = RELAX_get_segments(EEG);
        cmIdx = find(cm); un = cmIdx(~cmHasRef);
        U = struct('lo', pinfo.events.lo(un), 'hi', pinfo.events.hi(un), 'nChannels', pinfo.events.nChannels(un), ...
                   'signAgreement', nan(numel(un),1), 'amplitudeSpread', nan(numel(un),1), ...
                   'fractionNear', nan(numel(un),1), 'medianAmplitude', nan(numel(un),1));
        for j = 1:numel(un)
            seg = find(segments(:,1) <= U.lo(j) & segments(:,2) >= U.hi(j), 1);
            if isempty(seg); continue; end
            [~, u] = RELAX_step_uniformity(EEG.data, U.lo(j), U.hi(j), W, segments(seg,1), segments(seg,2));
            U.signAgreement(j) = u.signAgreement; U.amplitudeSpread(j) = u.amplitudeSpread;
            U.fractionNear(j) = u.fractionNear; U.medianAmplitude(j) = u.medianAmplitude;
        end
        r.unmatchedCommonMode = U;

        % 4. figures
        [~, base] = fileparts(opts.files{k});
        candidateFigures(EEG, r.ref{1}, opts.maxFigures, fullfile(opts.figDir, base), base);

        % 5. subtraction arm (only once Task 5 exists)
        if canSubtract
            r.subtract = subtractionArm(EEG, r.ref{1}, detectOpts, opts, slack);
            r.subtract.nCommonModeBefore = r.perChannel.nCommonModeEvents;
        end

        results{end+1} = r; %#ok<AGROW>
        save(outFile, 'results', 'opts', '-v7.3');
    end
    summaryFigure(results, fullfile(opts.figDir, 'summary_uniformity.png'));
    printSummary(results);
    save(outFile, 'results', 'opts', '-v7.3');
end

function dels = montageDeletions(montageFile)
    dels = {'E257'};
    if ~isempty(montageFile) && exist(montageFile, 'file') == 2
        idx = readmatrix(montageFile); idx = idx(~isnan(idx));
        drop = setdiff(1:257, idx(:)');
        dels = [arrayfun(@(n) sprintf('E%d', n), drop(:), 'UniformOutput', false); dels(:)];
    end
end

function printCandidates(c, srate)
    for k = 1:numel(c.sample)
        fprintf('  cand %2d: t = %8.2f s  amp %7.1f uV  rise %5.1f ms  agree %.2f  spread %.2f  near %.2f  z %.1f\n', ...
            k, c.sample(k) / srate, c.amplitude(k), c.riseMs(k), c.signAgreement(k), ...
            c.amplitudeSpread(k), c.fractionNear(k), c.z(k));
    end
end

function candidateFigures(EEG, R, maxFigures, prefix, name)
    c = R.candidates;
    if isempty(c.sample); return; end
    [~, order] = sort(abs(c.amplitude), 'descend');
    usedIdx = find(ismember({EEG.chanlocs.labels}, R.channelInfo.usedLabels));
    rng(1); show = usedIdx(randperm(numel(usedIdx), min(12, numel(usedIdx))));
    for k = order(1:min(maxFigures, numel(order)))'
        lo = max(1, c.onset(k) - EEG.srate); hi = min(EEG.pnts, c.offset(k) + EEG.srate);
        t = ((lo:hi) - c.sample(k)) / EEG.srate * 1000;
        fig = figure('Visible', 'off', 'Position', [0 0 900 700]);
        subplot(2,1,1); hold on
        for ch = show
            x = double(EEG.data(ch, lo:hi)); plot(t, x - median(x(1:EEG.srate)), 'Color', [0.7 0.7 0.7]);
        end
        rt = double(R.refTrace(lo:hi)); plot(t, rt - median(rt(1:EEG.srate)), 'k', 'LineWidth', 2);
        xline((c.onset(k) - c.sample(k)) / EEG.srate * 1000, 'r--'); xline((c.offset(k) - c.sample(k)) / EEG.srate * 1000, 'r--');
        xlabel('ms from mid-crossing'); ylabel('uV (demeaned over first 1 s)');
        title(sprintf('%s cand %d: amp %.0f uV, rise %.0f ms, agree %.2f, spread %.2f, near %.2f', ...
            strrep(name, '_', '\_'), k, c.amplitude(k), c.riseMs(k), c.signAgreement(k), c.amplitudeSpread(k), c.fractionNear(k)));
        subplot(2,1,2);
        a = R.perChannelAmplitudeAll(:, k); a = sort(a(isfinite(a)));
        bar(a); yline(median(a), 'r'); xlabel('channel (sorted)'); ylabel('level change, uV');
        print(fig, sprintf('%s_cand%02d.png', prefix, k), '-dpng', '-r90'); close(fig);
    end
end

function summaryFigure(results, path)
    fig = figure('Visible', 'off', 'Position', [0 0 700 500]); hold on
    for k = 1:numel(results)
        c = results{k}.ref{1}.candidates;
        scatter(c.amplitudeSpread, c.signAgreement, 30, 'b', 'filled');
        u = results{k}.unmatchedCommonMode;
        scatter(u.amplitudeSpread, u.signAgreement, 30, 'r');
    end
    xlim([0 2]); ylim([0.4 1.02]); xline(0.5, 'k--'); yline(0.9, 'k--');
    xlabel('amplitudeSpread (1.4826 MAD / |median|)'); ylabel('signAgreement');
    legend({'reference candidates', 'unmatched common-mode events'}, 'Location', 'southeast');
    title('Uniformity of level change across channels');
    print(fig, path, '-dpng', '-r90'); close(fig);
end

function S = subtractionArm(EEG, R, detectOpts, opts, slack)
    c = R.candidates; keep = R.gated;
    refSteps = structfun(@(f) f(keep), rmfield(c, 'accepted'), 'UniformOutput', false);
    [EEG2, sinfo] = RELAX_subtract_reference_pops(EEG, refSteps, double(R.refTrace));
    [~, after] = RELAX_detect_reference_pops(EEG2, 'MinSignAgreement', 0, 'MaxAmplitudeSpread', Inf, ...
        'Detect', [detectOpts, {'MinAmplitudeMicrovolts', opts.floors(1)}]);
    [~, pAfter] = RELAX_detect_steps(EEG2, detectOpts{:}, 'MinAmplitudeMicrovolts', 30, ...
        'CommonModeChannelFraction', opts.commonModeFraction);
    cm = pAfter.events.commonMode;
    [stillThere, ~] = RELAX_match_intervals(refSteps.onset, refSteps.offset, pAfter.events.lo(cm), pAfter.events.hi(cm), slack);
    S = struct('info', sinfo, 'nSubtracted', numel(refSteps.sample), ...
               'residualCandidates', after.candidates, 'nResidualGated', sum(after.candidates.signAgreement >= 0.9 & after.candidates.amplitudeSpread <= 0.5), ...
               'nCommonModeBefore', NaN, 'nCommonModeAfter', pAfter.nCommonModeEvents, ...
               'nStepsAfter', pAfter.nSteps, 'nSubtractedStillCommonMode', sum(stillThere));
    fprintf('Subtraction: %d reference pop(s) removed on %d channel(s); afterwards %d reference candidate(s) (%d gated), %d common-mode event(s) (%d at subtracted sites), %d per-channel steps\n', ...
        S.nSubtracted, sinfo.nChannelsCorrected, numel(after.candidates.sample), S.nResidualGated, ...
        S.nCommonModeAfter, S.nSubtractedStillCommonMode, S.nStepsAfter);
end

function printSummary(results)
    fprintf('\n===== SUMMARY (floor %g uV, gate agree >= 0.9 & spread <= 0.5)\n', results{1}.floors(1));
    fprintf('%-32s %6s %6s %6s %6s %6s %8s %8s\n', 'file', 'min', 'cand', 'gated', 'CMev', 'CM=ref', 'medAgr', 'medSpr');
    for k = 1:numel(results)
        r = results{k}; c = r.ref{1}.candidates; g = r.ref{1}.gated;
        fprintf('%-32s %6.1f %6d %6d %6d %6d %8.2f %8.2f\n', r.file, r.minutes, numel(c.sample), sum(g), ...
            r.perChannel.nCommonModeEvents, sum(r.match.commonModeHasRef), ...
            median(c.signAgreement(g)), median(c.amplitudeSpread(g)));
        if isfield(r, 'subtract')
            fprintf('%-32s   after subtraction: %d residual candidate(s), %d common-mode event(s)\n', '', ...
                numel(r.subtract.residualCandidates.sample), r.subtract.nCommonModeAfter);
        end
    end
end

function s = setDefault(s, f, v)
    if ~isfield(s, f) || isempty(s.(f)); s.(f) = v; end
end
