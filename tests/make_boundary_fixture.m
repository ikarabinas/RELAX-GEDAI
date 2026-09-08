function fixturePath = make_boundary_fixture()
% Build a small cropped real recording that contains boundary events, for the
% integration checks. Crops around the first boundary of a real input file so the
% fixture spans a dozen or so of the ~8 s inter-train segments.
%
% Returns the path to the written .set file.

    sourceDir  = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/';
    sourceFile = 'm121_dlpfc_day1_treatment_20221226_034250_interTrainIntervals.set';

    fixtureDir = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    if ~exist(fixtureDir, 'dir'); mkdir(fixtureDir); end
    fixturePath = fullfile(fixtureDir, 'boundary_fixture.set');

    if exist(fixturePath, 'file')
        fprintf('Fixture already exists at %s\n', fixturePath);
        return
    end

    EEG = pop_loadset('filename', sourceFile, 'filepath', sourceDir);

    segments = RELAX_get_segments(EEG);
    assert(size(segments,1) > 2, 'source file must contain boundary events');

    % Start 5 s before the first boundary and keep 120 s, which spans the first
    % boundary and roughly a dozen of the short inter-train segments after it.
    firstBoundary = segments(1,2);
    startSample = max(1, firstBoundary - 5*EEG.srate);
    endSample   = min(EEG.pnts, startSample + 120*EEG.srate - 1);

    EEG = pop_select(EEG, 'point', [startSample endSample]);
    EEG = eeg_checkset(EEG);

    croppedSegments = RELAX_get_segments(EEG);
    fprintf('Fixture: %d samples, %d segments, %d boundaries.\n', ...
        EEG.pnts, size(croppedSegments,1), size(croppedSegments,1)-1);
    assert(size(croppedSegments,1) > 5, 'cropped fixture must retain several segments');

    EEG = pop_saveset(EEG, 'filename', 'boundary_fixture.set', 'filepath', fixtureDir);
    fprintf('Wrote fixture to %s\n', fixturePath);
end
