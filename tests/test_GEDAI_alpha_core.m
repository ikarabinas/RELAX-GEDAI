function test_GEDAI_alpha_core()
% GEDAI_alpha_refcov builds the referenced leadfield covariance GEDAI is given, selecting rows by
% channel label. GEDAI_alpha_run makes one broadband-only GEDAI call and returns the off arm, the
% on arm and their difference. GEDAI_alpha_dipfit scores how well a topography is a single dipole.

    p0 = path;
    restorePath = onCleanup(@() path(p0));
    GEDAI_alpha_paths();

    rs = RandStream('mt19937ar', 'Seed', 3);
    nch = 32;
    nsrc = 200;
    L = randn(rs, nch, 3 * nsrc);
    labels = arrayfun(@(k) sprintf('E%d', k), 1:nch, 'UniformOutput', false);

    leadfieldDir = tempname;
    mkdir(leadfieldDir);
    removeDir = onCleanup(@() rmdir(leadfieldDir, 's'));
    ch_names = labels;
    save(fullfile(leadfieldDir, 'zz9_leadfield.mat'), 'L', 'ch_names');

    shuffled = labels(randperm(rs, nch));
    [refCOV, Lref] = GEDAI_alpha_refcov(shuffled, leadfieldDir, 'zz9');
    [~, rows] = ismember(shuffled, labels);
    Lr = L(rows, :);
    Lr = Lr - sum(Lr, 1) / (nch + 1);
    assert(isequal(size(Lref), [nch, 3 * nsrc]), '(a) Lref must be nch x 3*nsrc');
    assert(max(abs(Lref - Lr), [], 'all') < 1e-10, '(a) Lref must be the row-selected leadfield minus sum/(nch+1)');
    assert(max(abs(refCOV - Lr * Lr'), [], 'all') < 1e-10, '(a) refCOV must equal Lr*Lr'' in the order of the requested labels');

    try
        GEDAI_alpha_refcov([labels, {'E99'}], leadfieldDir, 'zz9');
        error('test:noError', '(b) a label missing from the leadfield must throw');
    catch err
        assert(strcmp(err.identifier, 'GEDAI_alpha:leadfieldChannels'), '(b) wrong error: %s', err.identifier);
    end

    gofDipole = GEDAI_alpha_dipfit(Lref(:, 31:33) * [1; -2; 0.5], Lref);
    assert(gofDipole > 0.999, '(c) a single free-orientation dipole must fit, gof = %g', gofDipole);
    gofNoise = GEDAI_alpha_dipfit(randn(rs, nch, 1), Lref);
    assert(gofNoise < 0.5, '(c) white-noise topography must not fit a dipole, gof = %g', gofNoise);

    opts = GEDAI_alpha_figures('opts');
    srate = 250;
    pnts = 120 * srate;
    [b, a] = butter(3, [1 40] / (srate / 2));
    rank = 8;
    U = randn(rs, nch, rank);
    Lbrain = U * randn(rs, rank, 3 * nsrc) + 0.01 * randn(rs, nch, 3 * nsrc);
    ch_names = labels;
    L = Lbrain;
    save(fullfile(leadfieldDir, 'zz8_leadfield.mat'), 'L', 'ch_names');
    sourceCols = [4 5 6; 100 101 102; 301 302 303];
    brain = zeros(nch, pnts);
    for s = 1:3
        activity = filter(b, a, randn(rs, 3, pnts), [], 2);
        brain = brain + Lbrain(:, sourceCols(s, :)) * activity;
    end
    brainRms = sqrt(mean(brain(:) .^ 2));
    Uref = U - sum(U, 1) / (nch + 1);
    pattern = randn(rs, nch, 1);
    pattern = pattern - Uref * (Uref \ pattern);
    pattern = pattern - sum(pattern) / (nch + 1);
    pattern = pattern / sqrt(mean(pattern .^ 2));
    patternTime = randn(rs, 1, pnts);
    data = brain + 5 * brainRms * pattern * patternTime;

    EEG = eeg_emptyset;
    EEG.data = data;
    EEG.nbchan = nch;
    EEG.pnts = pnts;
    EEG.trials = 1;
    EEG.srate = srate;
    EEG.xmin = 0;
    EEG.xmax = (pnts - 1) / srate;
    EEG.times = (0:pnts-1) / srate * 1000;
    EEG.chanlocs = struct('labels', labels);
    EEG.event = struct('type', 'boundary', 'latency', 60 * srate + 0.5, 'duration', 0);
    EEG.urevent = [];
    EEG = eeg_checkset(EEG, 'eventconsistency');

    refCOVrun = GEDAI_alpha_refcov(labels, leadfieldDir, 'zz8');
    A = GEDAI_alpha_run(EEG, refCOVrun, opts);

    assert(isequal(size(A.off), size(EEG.data)), '(d) off must have the input size');
    assert(isequal(size(A.on), size(A.off)) && isequal(size(A.removed), size(A.off)), '(d) on and removed must match off');
    scale = max(abs(A.off(:)));
    assert(max(abs(A.off - A.on - A.removed), [], 'all') < 1e-6 * scale, '(d) removed must equal off - on');
    assert(isscalar(A.threshold), '(d) threshold must be a scalar');
    assert(isscalar(A.enova) && isscalar(A.seconds), '(d) enova and seconds must be scalars');

    u = pattern / norm(pattern);
    offPower = mean((u' * A.off) .^ 2);
    removedPower = mean((u' * A.removed) .^ 2);
    assert(removedPower > 0.5 * offPower, ...
        '(e) positive control: removed must hold more than 50%% of the rank-1 pattern power, held %.3g of %.3g', ...
        removedPower, offPower);
end
