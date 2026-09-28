function test_RELAX_protect_eye_channels()
% Eye electrodes are exempt from pop-based rejection: a periocular electrode that trips a
% bad-contact criterion is far more likely to show eye movement the mask missed than to be
% a dead contact, and the pipeline's blink / HEOG stages cannot run without them (m407 lost
% all 14 of them and crashed). Its steps are left alone instead, like eye-movement steps.

    n = 200; rng(3);
    ch = randi(4, n, 1);
    steps = struct('channel', ch, 'label', {arrayfun(@(k) sprintf('E%d', k), ch, 'UniformOutput', false)}, ...
                   'sample', (1:n)' * 40, 'onset', (1:n)' * 40 - 5, 'offset', (1:n)' * 40 + 5, ...
                   'commonMode', false(n,1), 'eyeMovement', rand(n,1) < 0.2);
    info = struct('badChannels', {{'E3'; 'E1'; 'E4'}}, ...
                  'badChannelReason', {{'30.0 steps/min exceeds 6'; '12.0 steps/min exceeds 6'; '0.30 of samples in transition exceeds 0.02'}});
    eyeLabels = {'E1'; 'E10'; 'E37'};

    %% 1. an eye electrode leaves the rejection list, keeping the order of the rest
    [s1, i1] = RELAX_protect_eye_channels(steps, info, eyeLabels);
    assert(isequal(i1.badChannels, {'E3'; 'E4'}), 'E1 must be removed from the rejection list');
    assert(isequal(i1.badChannelReason, info.badChannelReason([1 3])), 'reasons must stay aligned with channels');
    assert(isequal(i1.protectedEyeChannels, {'E1'}) && iscell(i1.protectedEyeChannels));

    %% 2. its steps are all marked as left-alone (eyeMovement), others untouched
    onE1 = steps.channel == 1;
    assert(all(s1.eyeMovement(onE1)), 'every step on the protected electrode is left alone');
    assert(isequal(s1.eyeMovement(~onE1), steps.eyeMovement(~onE1)), 'other channels'' flags unchanged');
    assert(islogical(s1.eyeMovement) && isequal(size(s1.eyeMovement), [n 1]));

    %% 3. nothing to protect: identical output, 0x1 protected list
    info2 = struct('badChannels', {{'E3'}}, 'badChannelReason', {{'x'}});
    [s2, i2] = RELAX_protect_eye_channels(steps, info2, eyeLabels);
    assert(isequal(s2, steps) && isequal(i2.badChannels, {'E3'}));
    assert(iscell(i2.protectedEyeChannels) && isequal(size(i2.protectedEyeChannels), [0 1]));

    %% 4. empty rejection list and steps without an eyeMovement field
    steps3 = rmfield(steps, 'eyeMovement');
    info3 = struct('badChannels', {cell(0,1)}, 'badChannelReason', {cell(0,1)});
    [s3, i3] = RELAX_protect_eye_channels(steps3, info3, eyeLabels);
    assert(isequal(size(i3.badChannels), [0 1]) && isequal(size(i3.protectedEyeChannels), [0 1]));
    assert(isfield(s3, 'eyeMovement') && ~any(s3.eyeMovement) && isequal(size(s3.eyeMovement), [n 1]));

    %% 5. a protected electrode with no steps is still reported
    info4 = struct('badChannels', {{'E10'}}, 'badChannelReason', {{'0.30 of samples in transition exceeds 0.02'}});
    [~, i4] = RELAX_protect_eye_channels(steps, info4, eyeLabels);
    assert(isequal(i4.protectedEyeChannels, {'E10'}) && isempty(i4.badChannels));

    fprintf('test_RELAX_protect_eye_channels: PASS\n');
end
