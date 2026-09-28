function test_RELAX_flag_eye_steps()
% Eye-movement flagging: a multi-channel event whose eye electrodes (two of them, unless
% the event is common-mode) have their mid-crossing inside an eye-movement period, and
% none of whose members is over the amplitude ceiling, is skipped as a whole; anything
% else is left alone.

    function s = mk(channel, label, onset, offset, eventId)
        s = struct('channel', channel(:), 'label', {label(:)}, 'sample', onset(:) + 1, ...
                   'onset', onset(:), 'offset', offset(:), 'amplitude', 100 * ones(numel(onset),1), ...
                   'eventId', eventId(:), 'commonMode', false(numel(onset),1));
    end
    pnts = 10000;
    blink = {'E8', 'E9', 'E25'};
    mask = false(1, pnts); mask(4800:5400) = true;

    % 1. three-channel cluster containing E8 and E9, inside the mask -> all three flagged
    s = mk([3 8 9], {'E3','E8','E9'}, [5000 5005 4990], [5030 5040 5020], [1 1 1]);
    [out, info] = RELAX_flag_eye_steps(s, mask, blink);
    assert(islogical(out.eyeMovement) && isequal(size(out.eyeMovement), [3 1]));
    assert(all(out.eyeMovement), 'every member of an eye cluster is flagged');
    assert(info.nEyeSteps == 3 && info.nEyeEvents == 1);

    % 2. the same cluster outside the mask -> nothing
    out = RELAX_flag_eye_steps(s, false(1, pnts), blink);
    assert(~any(out.eyeMovement), 'outside eye periods nothing is flagged');

    % 3. a lone step on a blink electrode inside the mask is a pop, not a saccade
    s1 = mk(8, {'E8'}, 5000, 5030, 1);
    out = RELAX_flag_eye_steps(s1, mask, blink);
    assert(~out.eyeMovement, 'single-channel steps are never flagged');

    % 4. a cluster with no blink electrode inside the mask is left alone
    s2 = mk([3 4 5], {'E3','E4','E5'}, [5000 5005 4990], [5030 5040 5020], [1 1 1]);
    out = RELAX_flag_eye_steps(s2, mask, blink);
    assert(~any(out.eyeMovement));

    % 5. only the event whose mid-crossings sit in the mask is flagged; the last sample of
    %    the mask counts
    s3 = mk([8 9 8 9], {'E8','E9','E8','E9'}, [2000 2010 5390 5395], [2030 2040 5420 5430], [1 1 2 2]);
    s3.sample(3:4) = [5399; 5400];
    out = RELAX_flag_eye_steps(s3, mask, blink);
    assert(isequal(out.eyeMovement, [false; false; true; true]));

    % 6. empty input keeps 0x1 shapes; other fields untouched
    e = mk(zeros(0,1), cell(0,1), zeros(0,1), zeros(0,1), zeros(0,1));
    [out, info] = RELAX_flag_eye_steps(e, mask, blink);
    assert(isequal(size(out.eyeMovement), [0 1]) && islogical(out.eyeMovement));
    assert(info.nEyeSteps == 0 && info.nEyeEvents == 0);
    assert(isequal(out.sample, e.sample));

    % 7. one eye electrode is not enough for a local cluster (a 2-channel event where one
    %    channel happens to be E8 and the other is anywhere on the head is a pop pair),
    %    but a common-mode event only needs one
    s4 = mk([3 8], {'E3','E8'}, [5000 5005], [5030 5040], [1 1]);
    out = RELAX_flag_eye_steps(s4, mask, blink);
    assert(~any(out.eyeMovement), 'a local cluster with a single eye electrode is not flagged');
    s4.commonMode(:) = true;
    out = RELAX_flag_eye_steps(s4, mask, blink);
    assert(all(out.eyeMovement), 'a common-mode event with one eye electrode is flagged');

    % 8. the eye electrodes' mid-crossings must be inside the mask; overlap of the
    %    transition interval [onset, offset] alone does not count, and a non-eye member's
    %    mid-crossing does not stand in for them
    s5 = mk([3 8 9], {'E3','E8','E9'}, [5350 5390 5395], [5500 5500 5500], [1 1 1]);
    s5.sample = [5360; 5450; 5460];
    out = RELAX_flag_eye_steps(s5, mask, blink);
    assert(~any(out.eyeMovement), 'transition overlap without an eye mid-crossing in the mask is not eye movement');
    s5.sample(2) = 5395;
    out = RELAX_flag_eye_steps(s5, mask, blink);
    assert(all(out.eyeMovement), 'one eye mid-crossing in the mask is enough once two eye electrodes are present');

    % 9. a member over the amplitude ceiling (default 500 uV) makes the event a pop
    %    cluster, whichever channel it is on; the ceiling is configurable
    s6 = mk([3 8 9], {'E3','E8','E9'}, [5000 5005 4990], [5030 5040 5020], [1 1 1]);
    s6.amplitude(2) = -800;
    out = RELAX_flag_eye_steps(s6, mask, blink);
    assert(~any(out.eyeMovement), 'an 800 uV step on an eye electrode is not a saccade');
    s6.amplitude = [800; 100; 100];
    out = RELAX_flag_eye_steps(s6, mask, blink);
    assert(~any(out.eyeMovement), 'an 800 uV step on a non-eye member is not part of a saccade either');
    out = RELAX_flag_eye_steps(s6, mask, blink, 'MaxAmplitudeMicrovolts', 1000);
    assert(all(out.eyeMovement), 'the ceiling is configurable');
    s6.amplitude(:) = 500;
    out = RELAX_flag_eye_steps(s6, mask, blink);
    assert(all(out.eyeMovement), 'exactly at the ceiling is still eye movement');

    fprintf('test_RELAX_flag_eye_steps: PASS\n');
end
