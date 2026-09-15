% RELAX_channel_budget() - Running budget shared by the channel-rejection screens.
%
% Usage:
%   >> EEG          = RELAX_channel_budget(EEG, 'init', maxProportion);
%   >> [fits, EEG]  = RELAX_channel_budget(EEG, 'spend', stageName, worstFirstLabels);
%
% Inputs:
%   EEG              - EEGLAB EEG structure. 'init' reads .nbchan; 'spend' reads and
%                       writes EEG.RELAX.ChannelBudget.
%   maxProportion    - proportion of the ORIGINAL montage that all participating screens
%                       may delete between them. maxRejectable = floor(maxProportion*nbchan).
%   stageName        - char, a valid MATLAB field name ('pop', 'impedance', 'PREP'), used
%                       to attribute the spend in .spentBy.
%   worstFirstLabels - cell array of channel labels, ordered WORST FIRST. When the budget
%                       cannot cover the whole list it is truncated from the end, so the
%                       channels sacrificed are the least-bad ones. 'spend' truncates by
%                       keeping a PREFIX of this list, so the caller is responsible for
%                       the ordering; this function does not re-sort it.
%
% Outputs:
%   fits - the prefix of worstFirstLabels the budget allows. Always a cell array: a 1x0
%          empty cell when nothing fits (not the 0x0 empty cell a bare {} literal would
%          be), so callers can pass it straight to pop_select without a type test.
%   EEG  - EEG with .RELAX.ChannelBudget updated.
%
% Call 'init' immediately after EEG.allchan = EEG.chanlocs in the wrapper: that line
% records the reference montage the proportion is taken against, so the budget has to be
% created against the same snapshot, before any screen has removed a channel.
%
% IMPORTANT - worst-first is a caller obligation, and RELAX_impedance_bad_channels does
% NOT meet it: it selects the worst-N candidates internally but then sorts its return
% value back into ascending channel-index order (see its 'badIdx = sort(badIdx);' before
% building rejectLabels), so its output is channel-index order, not worst-first. Passing
% that output to 'spend' for truncation would therefore keep the lowest-numbered channels
% rather than the worst ones. The designed path is the other direction: pass the
% remaining budget (maxRejectable - nSpent) to RELAX_impedance_bad_channels as its 4th
% argument, maxRejectableCount, so it selects the worst N itself within the allowance;
% 'spend' truncation is then unreachable for that screen.
%
% Screens that do NOT participate, deliberately:
%   - the bridged-channel screen. A bridge is a physical short; a retained cluster member
%     holds a smeared mixture while its chanlocs entry claims a precise location, so
%     capping it would preserve channels known to be invalid.
%   - RELAX_excluding_channels_and_epoching, which keeps its own independent later cap.
%     Bringing it in means touching three sibling variants and is a separate change.

function [varargout] = RELAX_channel_budget(EEG, op, varargin)

    switch lower(op)
        case 'init'
            maxProportion = varargin{1};
            budget = struct( ...
                'nOriginalChannels', EEG.nbchan, ...
                'maxProportion',     maxProportion, ...
                'maxRejectable',     floor(maxProportion * EEG.nbchan), ...
                'nSpent',            0, ...
                'spentBy',           struct(), ...
                'exhausted',         false);
            EEG.RELAX.ChannelBudget = budget;  % leaves any other EEG.RELAX field alone
            varargout = {EEG};

        case 'spend'
            stageName        = varargin{1};
            worstFirstLabels = varargin{2};

            if ~isfield(EEG, 'RELAX') || ~isfield(EEG.RELAX, 'ChannelBudget')
                error('RELAX_channel_budget:notInitialised', ...
                    ['EEG.RELAX.ChannelBudget does not exist. Call ' ...
                     'RELAX_channel_budget(EEG, ''init'', maxProportion) before any screen spends.']);
            end

            if isempty(worstFirstLabels)
                worstFirstLabels = {};
            end
            worstFirstLabels = worstFirstLabels(:)';  % row, so the prefix slice is a row cell

            b = EEG.RELAX.ChannelBudget;
            remaining = max(0, b.maxRejectable - b.nSpent);
            nFits = min(numel(worstFirstLabels), remaining);
            fits = worstFirstLabels(1:nFits);

            b.nSpent = b.nSpent + nFits;
            b.spentBy.(stageName) = fits;
            if nFits < numel(worstFirstLabels)
                b.exhausted = true;
                warning('RELAX_channel_budget:exhausted', ...
                    ['Stage ''%s'' asked to reject %d channel(s) but only %d fit in the ' ...
                     'shared budget (%d of %d channels, %g of the montage). ' ...
                     'Rejecting the worst %d; the rest are recorded but retained.'], ...
                    stageName, numel(worstFirstLabels), nFits, ...
                    b.maxRejectable, b.nOriginalChannels, b.maxProportion, nFits);
            end
            EEG.RELAX.ChannelBudget = b;

            varargout = {fits, EEG};

        otherwise
            error('RELAX_channel_budget:unknownOp', ...
                'Unknown operation ''%s''. Expected ''init'' or ''spend''.', op);
    end
end
