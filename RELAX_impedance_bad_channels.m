% RELAX_impedance_bad_channels() - Select channels to reject on recorded electrode
%                                  impedance alone.
%
% Usage:
%   >> [rejectLabels, info] = RELAX_impedance_bad_channels(EEG.chanlocs, 500, 0.25);
%
% Inputs:
%   chanlocs      - EEG.chanlocs. The impedance is read from the .impedance field
%                   (kOhm), written at import by attach_mff_impedance in the
%                   TMS_response_variability repo. EEGLAB subsets chanlocs in
%                   lockstep with EEG.data, so the values still line up with the
%                   channels after earlier rejections.
%   maxKOhm       - Impedance limit. A channel is bad when its impedance is
%                   STRICTLY greater than this.
%   maxProportion - Safety cap: at most floor(maxProportion*nChan) channels may be
%                   rejected here. See below.
%   maxRejectableCount - (optional) explicit cap in channels, overriding
%                        floor(maxProportion*nChannels). Supplied by RELAX_channel_budget
%                        so this screen spends only the shared budget's remainder. Omit
%                        for the original single-screen behaviour.
%
% Outputs:
%   rejectLabels - Column cell array of channel labels to reject, {} if none.
%                  Labels rather than indices, so the caller can pass them
%                  straight to pop_select without tracking index shifts.
%   info         - Struct recording what the screen saw and did:
%                    .available       false when no channel carries a usable
%                                     impedance, in which case nothing is rejected
%                    .maxKOhm, .maxProportion   the limits applied
%                    .nChannels, .nWithImpedance, .nOverThreshold, .nRejected
%                    .capApplied      true when the cap held rejections back
%                    .rejectLabels, .rejectImpedances   the rejected channels
%                    .medianImpedance, .maxImpedance    over the known values
%
% This mirrors the impedance criterion in bad_channel_mask.m upstream, so the two
% screens agree on the same recording:
%
%   - the comparison is STRICTLY greater-than, so a channel exactly at the limit
%     is kept;
%   - an UNKNOWN (NaN) impedance is NOT bad. A recording whose .mff carried no
%     impedance block, or whose values failed to line up with the data channels,
%     must not have its whole montage rejected here;
%   - a recording with no usable impedance at all reports .available = false and
%     rejects nothing, loudly. Silently behaving like a clean recording would make
%     a missing impedance block indistinguishable from good electrodes.
%
% The cap exists because impedance is metadata: it is measured once, before the
% recording, and a cap-sized count of channels over the limit usually means the
% measurement itself is suspect rather than that most of the montage is dead.
% Rejecting the worst maxProportion and warning leaves the remaining channels to
% the data-driven screens downstream, which is a better failure mode than gutting
% the montage on one bad impedance sweep.

function [rejectLabels, info] = RELAX_impedance_bad_channels(chanlocs, maxKOhm, maxProportion, maxRejectableCount)

    % Optional 4th argument: an explicit cap in channels, from the shared budget in
    % RELAX_channel_budget. Overrides the internal floor(maxProportion*nChannels) so
    % this screen spends only what earlier screens left. The 3-argument form is
    % unchanged, which is what keeps the sibling wrappers byte-identical.
    if nargin < 4
        maxRejectableCount = [];
    end

    nChannels = numel(chanlocs);

    info = struct('available', false, ...
                  'maxKOhm', maxKOhm, ...
                  'maxProportion', maxProportion, ...
                  'nChannels', nChannels, ...
                  'nWithImpedance', 0, ...
                  'nOverThreshold', 0, ...
                  'nRejected', 0, ...
                  'capApplied', false, ...
                  'rejectLabels', {{}}, ...
                  'rejectImpedances', [], ...
                  'medianImpedance', NaN, ...
                  'maxImpedance', NaN);
    rejectLabels = {};

    if nChannels == 0
        return;
    end

    if ~isstruct(chanlocs) || ~isfield(chanlocs, 'impedance')
        warning('RELAX_impedance_bad_channels:noImpedanceField', ...
            ['chanlocs carries no .impedance field, so no channel can be screened on ' ...
             'impedance. Rejecting nothing. attach_mff_impedance writes this field at ' ...
             'import; a file imported without it has none.']);
        return;
    end

    % One value per channel, NaN where unknown or malformed. A chanlocs field can hold
    % anything, so anything that is not a real scalar becomes NaN rather than erroring:
    impedances = NaN(1, nChannels);
    for k = 1:nChannels
        thisValue = chanlocs(k).impedance;
        if isnumeric(thisValue) && isscalar(thisValue) && isreal(thisValue)
            impedances(k) = double(thisValue);
        end
    end

    isKnown = ~isnan(impedances);
    info.nWithImpedance = nnz(isKnown);

    if info.nWithImpedance == 0
        warning('RELAX_impedance_bad_channels:noImpedanceValues', ...
            ['No channel carries a known impedance (all %d are NaN), so no channel can ' ...
             'be screened on impedance. Rejecting nothing.'], nChannels);
        return;
    end

    info.available       = true;
    info.medianImpedance = median(impedances(isKnown));
    info.maxImpedance    = max(impedances(isKnown));

    % Strictly greater-than, and NaN > anything is false, so unknown impedances are
    % kept without needing a separate test:
    isBad = impedances > maxKOhm;
    info.nOverThreshold = nnz(isBad);

    if info.nOverThreshold == 0
        return;
    end

    badIdx = find(isBad);

    % Cap: keep the worst, drop the rest of the list. Sorting by impedance descending
    % makes "worst" mean the highest measured impedance:
    if isempty(maxRejectableCount)
        maxRejectable = floor(maxProportion * nChannels);
    else
        maxRejectable = maxRejectableCount;
    end
    if info.nOverThreshold > maxRejectable
        [~, worstFirst] = sort(impedances(badIdx), 'descend');
        badIdx = badIdx(worstFirst(1:maxRejectable));
        info.capApplied = true;
        warning('RELAX_impedance_bad_channels:tooManyOverThreshold', ...
            ['%d of %d channels are over %g kOhm, more than the %g cap (%d channels). ' ...
             'Rejecting only the %d worst and leaving the rest to the data-driven ' ...
             'screens - an impedance sweep this bad is more likely to be a bad ' ...
             'measurement than a mostly-dead montage. This file is worth inspecting.'], ...
            info.nOverThreshold, nChannels, maxKOhm, maxProportion, maxRejectable, maxRejectable);
    end

    badIdx = sort(badIdx);
    rejectLabels          = {chanlocs(badIdx).labels}';
    info.nRejected        = numel(badIdx);
    info.rejectLabels     = rejectLabels;
    info.rejectImpedances = impedances(badIdx);
end
