function [steps, info] = RELAX_protect_eye_channels(steps, info, eyeLabels)
% RELAX_protect_eye_channels - Exempt the eye electrodes from pop-based channel rejection.
%
%   [steps, info] = RELAX_protect_eye_channels(steps, info, eyeLabels)
%
% Removes every label in eyeLabels from info.badChannels / info.badChannelReason (order of
% the rest preserved) and marks all steps on those channels as eyeMovement, so
% RELAX_repair_steps leaves them alone. The removed labels are reported in
% info.protectedEyeChannels (0x1 cell when none).
%
% Why: a periocular electrode that exceeds the bad-contact criteria is far more likely to be
% showing eye movement the blink / HEOG mask missed than to be a dead contact (m407: HEOG
% mask 1.1 %, all 14 periocular electrodes "chronically popping"), and RELAX's own blink and
% HEOG stages cannot run without them. Leaving the steps unfilled also keeps the per-channel
% interpolation ceiling from being exceeded on such a channel.
%
% steps     - RELAX_detect_steps output (.label, optionally .eyeMovement)
% info      - struct with .badChannels / .badChannelReason (RELAX_rank_pop_channels output)
% eyeLabels - cell array of labels (BlinkElectrodes plus the HEOG pair the mask resolved)

    labels = steps.label(:);
    n = numel(labels);
    if ~isfield(steps, 'eyeMovement') || isempty(steps.eyeMovement)
        steps.eyeMovement = false(n, 1);
    end
    steps.eyeMovement = logical(steps.eyeMovement(:));

    bad = reshape(info.badChannels, [], 1);
    isEye = ismember(bad, eyeLabels(:));
    info.protectedEyeChannels = bad(isEye);
    info.badChannels = bad(~isEye);
    info.badChannelReason = reshape(info.badChannelReason, [], 1);
    info.badChannelReason = info.badChannelReason(~isEye);
    if isempty(info.protectedEyeChannels); info.protectedEyeChannels = cell(0, 1); end
    if isempty(info.badChannels); info.badChannels = cell(0, 1); info.badChannelReason = cell(0, 1); end

    steps.eyeMovement(ismember(labels, info.protectedEyeChannels)) = true;
end
