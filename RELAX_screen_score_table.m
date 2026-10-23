% RELAX_screen_score_table() - One row per file: the screening metrics, the screen verdict and the
%                              review-sample draw.
%
% Usage:
%   >> [screen, meta] = RELAX_screen_score_table(T, metrics, opts);
%
% Inputs:
%   T       - table from RELAX_QC_load_tables (File, Participant, Outcome)
%   metrics - [n x 1] struct array in the row order of T with the fields MuscleProportion,
%             MuscleStrength, BlinkRatio, NInterpolated, HoleDepth, HoleDeepestElectrode (string),
%             HoleElectrodes (string), Spacing, NInterpolatedLimit, NChannels and
%             MaxProportionDeleted (RELAX_screen_file_metrics), and InterpolatedReasons (string) with
%             NInterpolated_impedance, _pop, _bridged, _PREP, _extreme, _muscle and _unknown
%             (RELAX_screen_channel_reasons); NaN (strings "") for files that are not cleaned
%   opts    - .cutoff (a z-scored metric flags when z > cutoff), .holeDepthLimit (HoleDepth flags
%             when above it), .sampleSize (passing files drawn for review), .seed
%
% Outputs:
%   screen - one row per row of T, with the columns
%              File, Participant, PipelineOutcome;
%              for MuscleProportion and BlinkRatio: <M>, <M>_median, <M>_scale, <M>_z, <M>_flag;
%              MuscleStrength (a covariate: it has no z and never flags);
%              HoleDepth, HoleDepth_limit, HoleDepth_flag, HoleDeepestElectrode, HoleElectrodes;
%              NInterpolated, NInterpolated_limit (flags when the count reaches it),
%                NInterpolated_flag, NChannels, MaxProportionDeleted;
%              InterpolatedReasons (why each interpolated channel was removed, as 'muscle: E33 E34;
%                impedance: E5') and the counts NInterpolated_impedance, NInterpolated_pop,
%                NInterpolated_bridged, NInterpolated_PREP, NInterpolated_extreme,
%                NInterpolated_muscle and NInterpolated_unknown;
%              ScreenOutcome, FlaggedBy, InReviewSample, SheetStatus ("" until the sheets stage)
%                and SheetMessage
%   meta   - .cutoff, .holeDepthLimit, .sampleSize, .seed, .nCleaned, .created, .medianSpacing (of
%            the cleaned files, in the units of the electrode coordinates), .metricNames (the
%            decision and display order), .zMetricNames (the metrics with a robust z),
%            .zeroScaleMetrics (z metrics whose scale is not positive, so they cannot flag) and
%            .unscreeningMetrics (metrics whose NaN leaves a file unscreened)
%
% The median and scale come from the cleaned files only. The review sample uses its own
% seeded random stream, so the global rng is left alone.

function [screen, meta] = RELAX_screen_score_table(T, metrics, opts)

    zMetricNames = ["MuscleProportion", "BlinkRatio"];
    metricNames = [zMetricNames, "HoleDepth", "NInterpolated"];
    missingMakesUnscreened = [true true true false];
    nFiles = height(T);
    isCleaned = T.Outcome == "cleaned";

    screen = table(T.File, T.Participant, T.Outcome, 'VariableNames', {'File', 'Participant', 'PipelineOutcome'});
    values = NaN(nFiles, numel(metricNames));
    flagged = false(nFiles, numel(metricNames));
    zeroScale = false(1, numel(zMetricNames));
    for m = 1:numel(zMetricNames)
        name = zMetricNames(m);
        values(:, m) = [metrics.(name)]';
        [z, centre, scale] = RELAX_screen_robust_z(valuesOfCleanedOnly(values(:, m), isCleaned));
        zeroScale(m) = ~(scale > 0);
        flagged(:, m) = z > opts.cutoff & isCleaned;
        screen.(name) = values(:, m);
        screen.(name + "_median") = repmat(centre, nFiles, 1);
        screen.(name + "_scale") = repmat(scale, nFiles, 1);
        screen.(name + "_z") = z;
        screen.(name + "_flag") = flagged(:, m);
    end

    screen.MuscleStrength = [metrics.MuscleStrength]';

    depthColumn = numel(zMetricNames) + 1;
    values(:, depthColumn) = [metrics.HoleDepth]';
    flagged(:, depthColumn) = values(:, depthColumn) > opts.holeDepthLimit & isCleaned;
    screen.HoleDepth = values(:, depthColumn);
    screen.HoleDepth_limit = repmat(opts.holeDepthLimit, nFiles, 1);
    screen.HoleDepth_flag = flagged(:, depthColumn);
    screen.HoleDeepestElectrode = [metrics.HoleDeepestElectrode]';
    screen.HoleElectrodes = [metrics.HoleElectrodes]';

    countColumn = numel(zMetricNames) + 2;
    values(:, countColumn) = [metrics.NInterpolated]';
    limits = [metrics.NInterpolatedLimit]';
    flagged(:, countColumn) = values(:, countColumn) >= limits & isCleaned;
    screen.NInterpolated = values(:, countColumn);
    screen.NInterpolated_limit = limits;
    screen.NInterpolated_flag = flagged(:, countColumn);
    screen.NChannels = [metrics.NChannels]';
    screen.MaxProportionDeleted = [metrics.MaxProportionDeleted]';
    screen.InterpolatedReasons = [metrics.InterpolatedReasons]';
    for reasonName = ["impedance" "pop" "bridged" "PREP" "extreme" "muscle" "unknown"]
        screen.("NInterpolated_" + reasonName) = [metrics.("NInterpolated_" + reasonName)]';
    end

    [screen.ScreenOutcome, screen.FlaggedBy] = RELAX_screen_outcome(T.Outcome, values, flagged, ...
                                                                    metricNames, missingMakesUnscreened);
    screen.InReviewSample = drawReviewSample(screen.ScreenOutcome == "pass", opts.sampleSize, opts.seed);
    screen.SheetStatus = strings(nFiles, 1);
    screen.SheetMessage = strings(nFiles, 1);

    meta = struct('cutoff', opts.cutoff, 'holeDepthLimit', opts.holeDepthLimit, 'sampleSize', opts.sampleSize, ...
                  'seed', opts.seed, 'nCleaned', sum(isCleaned), 'created', datetime('now'), ...
                  'medianSpacing', median([metrics(isCleaned).Spacing], 'omitnan'), ...
                  'metricNames', metricNames, 'zMetricNames', zMetricNames, ...
                  'zeroScaleMetrics', zMetricNames(zeroScale), 'unscreeningMetrics', metricNames(missingMakesUnscreened));
end

function x = valuesOfCleanedOnly(x, isCleaned)
    x(~isCleaned) = NaN;
end

function inSample = drawReviewSample(isPass, sampleSize, seed)
    stream = RandStream('mt19937ar', 'Seed', seed);
    passing = find(isPass);
    chosen = passing(randperm(stream, numel(passing), min(sampleSize, numel(passing))));
    inSample = false(size(isPass));
    inSample(chosen) = true;
end
