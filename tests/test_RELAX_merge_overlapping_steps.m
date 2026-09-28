function test_RELAX_merge_overlapping_steps()
% Steps whose [onset, offset] intervals overlap or coincide become one step: the union
% interval, the amplitude re-measured on the trace with the detector's windows, the
% larger z, riseMs from the union. Non-overlapping steps pass through unchanged, in onset
% order, and every input step knows which output row it went to.

    srate = 1000; n = 6000; W = 150;
    t = 1:n;
    x = 60 * (t >= 1050) + 60 * (t >= 1150) ...    % a staircase bracketed twice
      + 70 * (t >= 3001) ...                       % a single clean step
      + 40 * (t >= 4992);                          % a step 8 samples from the segment end
    segments = [1 5000; 5001 n];
    ms = @(on, off) (off - on - 1) / srate * 1000;

    in = struct('sample', [1150; 1100; 3000; 3000; 4991], ...
                'onset',  [1140; 1040; 2999; 2999; 4990], ...
                'offset', [1160; 1160; 3001; 3001; 4992], ...
                'amplitude', [60; 120; 70; 70; 40], ...
                'riseMs', [ms(1140,1160); ms(1040,1160); ms(2999,3001); ms(2999,3001); ms(4990,4992)], ...
                'z', [5; 9; 6; 8; 4], 'channelZ', [7; 3; 6; 6; 5], 'extra', [1; 2; 3; 4; 5]);

    [out, groupOf] = RELAX_merge_overlapping_steps(in, x, segments, W, srate);

    %% 1. one row per run, in onset order; membership recorded
    assert(isequal(groupOf, [1; 1; 2; 2; 3]), 'groupOf %s', mat2str(groupOf'));
    for f = fieldnames(in)'
        assert(isequal(size(out.(f{1})), [3 1]), '%s must be 3 x 1', f{1});
    end
    assert(isequal(out.onset, [1040; 2999; 4990]) && isequal(out.offset, [1160; 3001; 4992]));

    %% 2. the staircase: union interval, amplitude 120 (not 60 + 120), z = max, lead row's extras
    assert(abs(out.amplitude(1) - 120) < 1e-9, 'amplitude %.1f', out.amplitude(1));
    assert(out.z(1) == 9 && out.channelZ(1) == 7);
    assert(abs(out.riseMs(1) - ms(1040, 1160)) < 1e-9);
    assert(out.sample(1) == 1100 && out.extra(1) == 2, 'other columns come from the member with the largest |z|');

    %% 3. exact duplicates collapse to one row that equals the lead member
    assert(out.amplitude(2) == 70 && out.z(2) == 8 && out.extra(2) == 4 && out.riseMs(2) == ms(2999,3001));

    %% 4. a lone step passes through untouched, and the post window is clipped at the segment end
    assert(isequal(structfun(@(f) f(3), out), structfun(@(f) f(5), in)));
    dup = in; dup.onset(5) = 4989; dup.offset(5) = 4992;              % overlap the last step with a wider copy
    dup = structfun(@(f) [f; f(5)], dup, 'UniformOutput', false); dup.onset(6) = 4990; dup.z(6) = 1;
    out4 = RELAX_merge_overlapping_steps(dup, x, segments, W, srate);
    assert(abs(out4.amplitude(3) - 40) < 1e-9, 'amplitude measured on x(4992:5000) - x(..:4989): %.1f', out4.amplitude(3));

    %% 5. empty in, empty (0 x 1) out
    e = structfun(@(f) zeros(0, 1), in, 'UniformOutput', false);
    [outE, groupE] = RELAX_merge_overlapping_steps(e, x, segments, W, srate);
    assert(isequal(size(outE.onset), [0 1]) && isequal(size(groupE), [0 1]));

    %% 6. touching intervals (offset == next onset) share a sample and merge; adjacent ones do not
    t = struct('sample', [100; 200], 'onset', [90; 110], 'offset', [110; 210], 'amplitude', [1; 1], ...
               'riseMs', [0; 0], 'z', [1; 1], 'channelZ', [1; 1]);
    [~, gT] = RELAX_merge_overlapping_steps(t, x, segments, W, srate);
    assert(isequal(gT, [1; 1]));
    t.onset(2) = 111;
    [~, gA] = RELAX_merge_overlapping_steps(t, x, segments, W, srate);
    assert(isequal(gA, [1; 2]));

    fprintf('test_RELAX_merge_overlapping_steps: PASS\n');
end
