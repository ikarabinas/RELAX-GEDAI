function test_RELAX_IC_muscle_slopes()
% Muscle-slope fit of ICA activations, whole-record and segment by segment.

    rng(11);
    srate = 500;
    segmentSeconds = 8;
    nSegments = 6;
    samplesPerSegment = segmentSeconds * srate;
    muscleLike = diff(randn(2, nSegments * samplesPerSegment + 1), 1, 2);
    brainLike = cumsum(randn(2, nSegments * samplesPerSegment), 2);
    brainLike = brainLike - movmean(brainLike, srate, 2);

    slopes = RELAX_IC_muscle_slopes([muscleLike; brainLike], srate, 50);
    assert(all(slopes(1:2) > 1) && all(slopes(3:4) < -1), ...
        'differenced noise should fit a rising slope and integrated noise a falling one');

    singleSegment = [1, nSegments * samplesPerSegment];
    assert(isequal(RELAX_IC_muscle_slopes(muscleLike, srate, 50, singleSegment), ...
        RELAX_IC_muscle_slopes(muscleLike, srate, 50)), ...
        'one segment spanning the record should equal the whole-record fit');

    starts = (0:nSegments - 1)' * samplesPerSegment + 1;
    segments = [starts, starts + samplesPerSegment - 1];
    offsets = repelem(30 * randn(2, nSegments), 1, samplesPerSegment);
    spliced = muscleLike(:, 1:nSegments * samplesPerSegment) + offsets;
    wholeRecord = RELAX_IC_muscle_slopes(spliced, srate, 50);
    bySegment = RELAX_IC_muscle_slopes(spliced, srate, 50, segments);
    assert(all(wholeRecord < -0.31), ...
        'steps between segments should swamp the 12-70 Hz muscle power in a whole-record fit');
    assert(all(bySegment > 1), ...
        'fitting within each segment should recover the muscle slope despite the steps');

    sliver = [segments; segments(end, 2) + 1, segments(end, 2) + srate];
    withSliver = RELAX_IC_muscle_slopes([spliced, 1000 * randn(2, srate)], srate, 50, sliver);
    assert(max(abs(withSliver - bySegment)) < 1e-12, ...
        'a segment under two seconds should be left out of the fit');

    threw = false;
    try
        RELAX_IC_muscle_slopes(spliced, srate, 50, [1, srate]);
    catch
        threw = true;
    end
    assert(threw, 'segments all too short for a spectrum should error, not return a slope');
end
