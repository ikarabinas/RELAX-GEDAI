function cache = make_qc_cache_fixture(nOk)
% Synthetic stage-1 cache shaped like RELAX_QC_extract's output: nOk aligned files, then
% one alignment-skipped file, one excluded recording and one never-processed file.
    if nargin < 1, nOk = 4; end
    rng(2);
    chanlocs = make_qc_chanlocs(32);
    freqs = 0:0.5:80;
    blinkTimes = (-250:250) / 500;
    n = nOk + 3;
    files = "f" + (1:n)';
    outcome = [repmat("cleaned", nOk + 1, 1); "excluded"; "notProcessed"];
    message = [repmat("", nOk + 1, 1); "excessive bridging"; ""];
    T = table(files, repmat("p1", n, 1), (1:n)', outcome, message, ...
              'VariableNames', {'File', 'Participant', 'FileNumber', 'Outcome', 'Message'});
    T.NChanImpedance = randi(3, n, 1); T.NChanPop = randi(2, n, 1); T.NChanBridge = nan(n, 1);
    T.NChanPREP = zeros(n, 1); T.NChanMuscle = zeros(n, 1); T.NChanExtreme = randi(2, n, 1);
    T.ExtremeExcisedSeconds = 30 * rand(n, 1); T.PopExcisedSeconds = zeros(n, 1);
    T.PopFractionSamplesInterpolated = 1e-3 * rand(n, 1); T.PopNSteps = randi(800, n, 1);
    T.BlinkElectrodeFallback = zeros(n, 1); T.WicaArtifactICsReduced = 0.4 * rand(n, 1);
    T.TableRawBlinkRatio = repmat({NaN}, n, 1); T.TableCleanedBlinkRatio = repmat({NaN}, n, 1);
    T.TableCleanedMuscleProportion = nan(n, 1); T.TableCleanedMuscleStrength = nan(n, 1);
    T{nOk + 2:n, {'NChanImpedance', 'NChanPop', 'NChanExtreme', 'ExtremeExcisedSeconds', 'PopNSteps', 'WicaArtifactICsReduced'}} = NaN;

    summaries = repmat(RELAX_QC_empty_summary(), n, 1);
    alpha = exp(-((freqs - 10) / 2).^2);
    eye = 5 * exp(-freqs / 3);
    muscle = 0.02 * (freqs / 80);
    shape = 100 * exp(-(blinkTimes / 0.1).^2);
    frontness = max([chanlocs.X]', 0);
    for i = 1:nOk
        s = summaries(i);
        s.status = "ok"; s.durationSeconds = 1000; s.nChannels = 32; s.chanlocs = chanlocs; s.inputSrate = 1000;
        s.freqs = freqs;
        s.psdBaseline = 1 ./ max(freqs, 0.5) + alpha + eye + muscle;
        s.psdRemoved = 0.9 * (eye + muscle);
        s.psdCleaned = s.psdBaseline - s.psdRemoved;
        s.psdInput = 1.5 * s.psdBaseline;
        s.varianceRemovedFraction = 0.3 * rand;
        s.blinkN = 20; s.blinkTimes = blinkTimes;
        s.blinkErpBaseline = shape; s.blinkErpCleaned = 0.1 * shape; s.blinkErpRemoved = 0.9 * shape;
        s.topoBaseline = 100 * frontness; s.topoCleaned = 0.1 * s.topoBaseline + randn(32, 1); s.topoRemoved = s.topoBaseline - s.topoCleaned;
        s.metricsBaseline.blinkRatio = 3 + rand; s.metricsCleaned.blinkRatio = 1 + 0.2 * rand;
        s.metricsBaseline.muscleProportion = 0.3 * rand; s.metricsCleaned.muscleProportion = 0.1 * rand;
        s.metricsBaseline.muscleStrength = 0.2 * rand; s.metricsCleaned.muscleStrength = 0.05 * rand;
        s.scores.blinkAmplitude = 100 * rand; s.scores.musclePower = rand; s.scores.popSteps = T.PopNSteps(i);
        summaries(i) = s;
    end
    if nOk >= 1, summaries(nOk).blinkN = 2; end
    summaries(nOk + 1).status = "skipped";
    summaries(nOk + 1).skipReason = "pnts differ: baseline 10, cleaned 9";

    examples = struct('file', {}, 'role', {}, 'windows', {});
    picks = struct('file', {}, 'role', {}, 'centres', {});
    roles = ["median", "worst blink"];
    for k = 1:min(2, nOk)
        windows = [exampleWindow("blink", chanlocs), exampleWindow("muscle", chanlocs), exampleWindow("pop", chanlocs)];
        if k == 2
            windows(3) = RELAX_QC_empty_window("pop");
            windows(3).note = "No repaired pop in this file";
        end
        examples(end+1) = struct('file', files(k), 'role', roles(k), 'windows', windows); %#ok<AGROW>
        picks(end+1) = struct('file', files(k), 'role', roles(k), 'centres', struct('blink', NaN, 'muscle', NaN, 'pop', NaN)); %#ok<AGROW>
    end

    cache.meta = struct('created', datetime('now'), 'myPath', "", 'opts', struct(), 'muscleSlopeThreshold', -0.31, ...
                        'maxChannelProportion', 0.25, 'minBlinksForGroup', 5);
    cache.tables = T;
    cache.summaries = summaries;
    cache.picks = picks;
    cache.examples = examples;
end

function w = exampleWindow(kind, chanlocs)
    w = RELAX_QC_empty_window(kind);
    w.available = true;
    w.channels = {chanlocs(1:6).labels};
    w.t = 100 + (0:2999) / 500;
    bump = 150 * exp(-((w.t - 103) / 0.15).^2);
    w.before = 10 * randn(6, 3000) + [bump; bump; 0.5 * bump; zeros(3, 3000)];
    w.after = w.before - [0.95 * bump; 0.95 * bump; 0.45 * bump; zeros(3, 3000)];
    w.removed = w.before - w.after;
    w.boundaries = 104.5;
    if kind == "pop"
        w.shade = [102.9 103.1];
        w.excised = [101.2 101.8; 104.6 105.4];
        w.note = "Grey: input band-passed without pop repair (reconstruction, not saved pipeline data)";
    end
end
