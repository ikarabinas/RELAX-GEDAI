% RELAX_screen_score_table() - One row per file: the four screening metrics, their robust
%                              z-scores, the screen verdict and the review-sample draw.
%
% Usage:
%   >> [screen, meta] = RELAX_screen_score_table(T, metrics, opts);
%
% Inputs:
%   T       - table from RELAX_QC_load_tables (File, Participant, Outcome)
%   metrics - [n x 1] struct array in the row order of T with the fields MuscleProportion,
%             MuscleStrength, BlinkRatio, NInterpolated; NaN for files that are not cleaned
%   opts    - .cutoff (a metric flags when z > cutoff), .sampleSize (passing files drawn for
%             review), .seed
%
% Outputs:
%   screen - one row per row of T: File, Participant, PipelineOutcome; for each metric <M>,
%            <M>_median, <M>_scale, <M>_z, <M>_flag; ScreenOutcome, FlaggedBy, InReviewSample,
%            SheetStatus ("" until the sheets stage) and SheetMessage
%   meta   - .cutoff, .sampleSize, .seed, .nCleaned, .created, .metricNames and
%            .zeroScaleMetrics (metrics whose scale is not positive, so they cannot flag) and
%            .unscreeningMetrics (metrics whose NaN leaves a file unscreened)
%
% The median and scale come from the cleaned files only. The review sample uses its own
% seeded random stream, so the global rng is left alone.

function [screen, meta] = RELAX_screen_score_table(T, metrics, opts)

    metricNames = ["MuscleProportion", "MuscleStrength", "BlinkRatio", "NInterpolated"];
    missingMakesUnscreened = [true false true true];
    nFiles = height(T);
    isCleaned = T.Outcome == "cleaned";

    screen = table(T.File, T.Participant, T.Outcome, 'VariableNames', {'File', 'Participant', 'PipelineOutcome'});
    values = NaN(nFiles, numel(metricNames));
    zScores = NaN(nFiles, numel(metricNames));
    zeroScale = false(1, numel(metricNames));
    for m = 1:numel(metricNames)
        name = metricNames(m);
        values(:, m) = [metrics.(name)]';
        [zScores(:, m), centre, scale] = RELAX_screen_robust_z(valuesOfCleanedOnly(values(:, m), isCleaned));
        zeroScale(m) = ~(scale > 0);
        screen.(name) = values(:, m);
        screen.(name + "_median") = repmat(centre, nFiles, 1);
        screen.(name + "_scale") = repmat(scale, nFiles, 1);
        screen.(name + "_z") = zScores(:, m);
        screen.(name + "_flag") = zScores(:, m) > opts.cutoff & isCleaned;
    end

    [screen.ScreenOutcome, screen.FlaggedBy] = RELAX_screen_outcome(T.Outcome, values, zScores, opts.cutoff, ...
                                                                    metricNames, missingMakesUnscreened);
    screen.InReviewSample = drawReviewSample(screen.ScreenOutcome == "pass", opts.sampleSize, opts.seed);
    screen.SheetStatus = strings(nFiles, 1);
    screen.SheetMessage = strings(nFiles, 1);

    meta = struct('cutoff', opts.cutoff, 'sampleSize', opts.sampleSize, 'seed', opts.seed, ...
                  'nCleaned', sum(isCleaned), 'created', datetime('now'), ...
                  'metricNames', metricNames, 'zeroScaleMetrics', metricNames(zeroScale), ...
                  'unscreeningMetrics', metricNames(missingMakesUnscreened));
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
