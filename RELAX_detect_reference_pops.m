% RELAX_detect_reference_pops() - Find pops on the recording reference electrode as steps
%                                 on the cross-channel median, and check that each one
%                                 really is the same level change on every channel.
%
% Usage:
%   >> [refSteps, refInfo] = RELAX_detect_reference_pops(EEG);
%   >> [refSteps, refInfo] = RELAX_detect_reference_pops(EEG, 'Detect', ...
%          {'MinAmplitudeMicrovolts', 30, 'WindowMs', 150}, 'MinSignAgreement', 0.9);
%
% Every channel is x_i(t) = s_i(t) - r(t), so a pop on the reference contact r is the same
% waveform on every channel and survives the cross-channel median
% (RELAX_reference_channel), while contact pops, movement and eye activity do not. The
% existing two-window detector (RELAX_detect_steps) is run on that one virtual channel;
% each candidate is then measured on the multichannel data (RELAX_step_uniformity) and
% kept only if the per-channel level changes agree in sign and amplitude. A same-sign
% event with widely varying amplitude (a ground-electrode pop, a partial movement) is a
% candidate on the median but fails the gate and is left to the per-channel stage.
% Accepted candidates whose intervals overlap (a staircase bracketed twice, the same
% bracket at two z values) are merged into one step (RELAX_merge_overlapping_steps), so
% that RELAX_subtract_reference_pops never corrects a transition twice.
%
% Run BEFORE any filtering, at the native rate, and before the per-channel pop detector.
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   Name-value pairs (all optional):
%     ExcludeLabels      - channels left out of the median and the uniformity vector
%     MinChannels        - fewer usable channels -> unavailable (default 8)
%     Detect             - cell of name-value pairs forwarded to RELAX_detect_steps
%                          (default {}: that function's own defaults)
%     MinSignAgreement   - gate on summary.signAgreement (default 0.9)
%     MaxAmplitudeSpread - gate on summary.amplitudeSpread (default 0.5)
%
% Outputs:
%   refSteps - struct of parallel N x 1 columns (0x1 when empty), on EEG's timeline, one
%              row per merged run of accepted candidates: .sample .onset .offset
%              .amplitude (post minus pre on the median, i.e. the change every channel
%              shows) .riseMs .z .channelZ .signAgreement .amplitudeSpread .fractionNear
%   refInfo  - .available .status .params .refChannel (the virtual 1-channel EEG)
%              .channelInfo (RELAX_reference_channel info) .detectInfo (RELAX_detect_steps
%              info) .minutes .nCandidates .nSteps (after merging) .nRejectedUniformity
%              .nMerged (accepted candidates absorbed into another) .stepsPerMinute
%              .perChannelAmplitude (nbchan x nSteps; NaN on channels not used)
%              .candidates (refSteps' columns for EVERY candidate plus .accepted,
%              .mergedInto (the refSteps row it went to; 0 when rejected), and
%              .perChannelAmplitude (nbchan x nCandidates; NaN on channels not used) -
%              the per-candidate counterpart of refInfo.perChannelAmplitude, needed
%              because merging can leave nSteps < nCandidates even when every
%              candidate is accepted (e.g. survey code running with the uniformity
%              gate off))

function [refSteps, refInfo] = RELAX_detect_reference_pops(EEG, varargin)
    p = inputParser;
    p.addParameter('ExcludeLabels', {});
    p.addParameter('MinChannels', 8);
    p.addParameter('Detect', {});
    p.addParameter('MinSignAgreement', 0.9);
    p.addParameter('MaxAmplitudeSpread', 0.5);
    p.parse(varargin{:});
    opt = p.Results;

    nb = EEG.nbchan;
    refSteps = emptyRefSteps();
    [ref, chanInfo] = RELAX_reference_channel(EEG, 'ExcludeLabels', opt.ExcludeLabels, ...
                                              'MinChannels', opt.MinChannels);
    refInfo = struct('available', false, 'status', chanInfo.status, 'params', opt, ...
                     'refChannel', ref, 'channelInfo', chanInfo, 'detectInfo', [], ...
                     'minutes', EEG.pnts / EEG.srate / 60, 'nCandidates', 0, 'nSteps', 0, ...
                     'nRejectedUniformity', 0, 'nMerged', 0, 'stepsPerMinute', 0, ...
                     'perChannelAmplitude', nan(nb, 0), 'candidates', emptyRefSteps());
    refInfo.candidates.accepted = false(0, 1);
    refInfo.candidates.mergedInto = zeros(0, 1);
    if ~chanInfo.available; return; end

    [s, d] = RELAX_detect_steps(ref, opt.Detect{:});
    refInfo.detectInfo = d;
    refInfo.status = d.status;
    if ~d.available; return; end
    refInfo.available = true;

    W = d.params.derived.W;
    segments = chanInfo.segments;
    usedIdx = find(ismember(reshape({EEG.chanlocs.labels}, [], 1), chanInfo.usedLabels));
    nCand = numel(s.sample);
    amp = nan(nb, nCand);
    agreement = nan(nCand, 1); spread = nan(nCand, 1); near = nan(nCand, 1);
    for k = 1:nCand
        [amp(usedIdx, k), agreement(k), spread(k), near(k)] = uniformity(s.onset(k), s.offset(k));
    end
    accepted = agreement >= opt.MinSignAgreement & spread <= opt.MaxAmplitudeSpread;

    cand = struct('sample', s.sample, 'onset', s.onset, 'offset', s.offset, ...
                  'amplitude', s.amplitude, 'riseMs', s.riseMs, 'z', s.z, 'channelZ', s.channelZ, ...
                  'signAgreement', agreement, 'amplitudeSpread', spread, 'fractionNear', near);
    acceptedSteps = structfun(@(f) f(accepted), cand, 'UniformOutput', false);
    [refSteps, groupOf] = RELAX_merge_overlapping_steps(acceptedSteps, ref.data, segments, W, EEG.srate);
    nSteps = numel(refSteps.onset);
    ampAccepted = amp(:, accepted);
    ampSteps = nan(nb, nSteps);
    for g = 1:nSteps
        members = find(groupOf == g);
        if isscalar(members)
            ampSteps(:, g) = ampAccepted(:, members);
        else
            [ampSteps(usedIdx, g), refSteps.signAgreement(g), refSteps.amplitudeSpread(g), ...
                refSteps.fractionNear(g)] = uniformity(refSteps.onset(g), refSteps.offset(g));
        end
    end
    cand.accepted = accepted;
    cand.mergedInto = zeros(nCand, 1);
    cand.mergedInto(accepted) = groupOf;
    cand.perChannelAmplitude = amp;

    refInfo.candidates = cand;
    refInfo.nCandidates = nCand;
    refInfo.nSteps = nSteps;
    refInfo.nRejectedUniformity = sum(~accepted);
    refInfo.nMerged = sum(accepted) - nSteps;
    refInfo.stepsPerMinute = nSteps / refInfo.minutes;
    refInfo.perChannelAmplitude = ampSteps;

    function [a, agree, spr, nr] = uniformity(onset, offset)
        seg = find(segments(:,1) <= onset & segments(:,2) >= offset, 1);
        [a, u] = RELAX_step_uniformity(EEG.data(usedIdx, :), onset, offset, W, ...
                                       segments(seg,1), segments(seg,2));
        agree = u.signAgreement; spr = u.amplitudeSpread; nr = u.fractionNear;
    end
end

function s = emptyRefSteps()
    z = zeros(0, 1);
    s = struct('sample', z, 'onset', z, 'offset', z, 'amplitude', z, 'riseMs', z, ...
               'z', z, 'channelZ', z, ...
               'signAgreement', z, 'amplitudeSpread', z, 'fractionNear', z);
end
