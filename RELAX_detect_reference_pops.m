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
%   refSteps - struct of parallel N x 1 columns (0x1 when empty), on EEG's timeline:
%              .sample .onset .offset .amplitude (post minus pre on the median, i.e. the
%              change every channel shows) .riseMs .z .channelZ .signAgreement
%              .amplitudeSpread .fractionNear
%   refInfo  - .available .status .params .refChannel (the virtual 1-channel EEG)
%              .channelInfo (RELAX_reference_channel info) .detectInfo (RELAX_detect_steps
%              info) .minutes .nCandidates .nSteps .nRejectedUniformity .stepsPerMinute
%              .perChannelAmplitude (nbchan x nSteps; NaN on channels not used)
%              .candidates (refSteps' columns for EVERY candidate plus .accepted)

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
                     'nRejectedUniformity', 0, 'stepsPerMinute', 0, ...
                     'perChannelAmplitude', nan(nb, 0), 'candidates', emptyRefSteps());
    refInfo.candidates.accepted = false(0, 1);
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
        seg = find(segments(:,1) <= s.onset(k) & segments(:,2) >= s.offset(k), 1);
        [a, u] = RELAX_step_uniformity(EEG.data(usedIdx, :), s.onset(k), s.offset(k), W, ...
                                       segments(seg,1), segments(seg,2));
        amp(usedIdx, k) = a;
        agreement(k) = u.signAgreement; spread(k) = u.amplitudeSpread; near(k) = u.fractionNear;
    end
    accepted = agreement >= opt.MinSignAgreement & spread <= opt.MaxAmplitudeSpread;

    cand = struct('sample', s.sample, 'onset', s.onset, 'offset', s.offset, ...
                  'amplitude', s.amplitude, 'riseMs', s.riseMs, 'z', s.z, 'channelZ', s.channelZ, ...
                  'signAgreement', agreement, 'amplitudeSpread', spread, 'fractionNear', near);
    refSteps = structfun(@(f) f(accepted), cand, 'UniformOutput', false);
    cand.accepted = accepted;

    refInfo.candidates = cand;
    refInfo.nCandidates = nCand;
    refInfo.nSteps = sum(accepted);
    refInfo.nRejectedUniformity = sum(~accepted);
    refInfo.stepsPerMinute = refInfo.nSteps / refInfo.minutes;
    refInfo.perChannelAmplitude = amp(:, accepted);
end

function s = emptyRefSteps()
    z = zeros(0, 1);
    s = struct('sample', z, 'onset', z, 'offset', z, 'amplitude', z, 'riseMs', z, ...
               'z', z, 'channelZ', z, 'signAgreement', z, 'amplitudeSpread', z, 'fractionNear', z);
end
