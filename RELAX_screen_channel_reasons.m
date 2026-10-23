% RELAX_screen_channel_reasons() - Why each interpolated electrode of a cleaned file was removed.
%
% Usage:
%   >> reasons = RELAX_screen_channel_reasons(EEG);
%
% Inputs:
%   EEG - a cleaned file's header or the full EEG. Needs allchan and RELAX.ListOfChannelsAfterRejections;
%         uses RELAXProcessingExtremeRejections and RELAX.PopScreen when present. An electrode listed by
%         more than one screen takes the first reason in the order impedance, pop, bridged, PREP, extreme,
%         muscle. An electrode no screen lists, or every electrode when the screens are absent, is "unknown".
%
% Outputs:
%   reasons.labels  - [1 x n] string, the interpolated electrodes (allchan labels not in
%                     RELAX.ListOfChannelsAfterRejections), in allchan order
%   reasons.reasons - [1 x n] string, the reason for each label. A pop reason names its sub-reason from
%                     RELAX.PopScreen, for example "pop (rate)", and is plain "pop" without one
%   reasons.summary - string, one "reason: labels" group per reason joined with "; " (labels joined with
%                     " "), in reason order. "" when nothing is interpolated
%   reasons.counts  - struct with fields impedance, pop, bridged, PREP, extreme, muscle and unknown, the
%                     number of interpolated electrodes per reason; every pop sub-reason counts as pop
%   reasons.groups  - struct array with fields reason (string) and labels ([1 x k] string), one per
%                     non-empty full reason, in reason order and, within pop, in order of the first label

function reasons = RELAX_screen_channel_reasons(EEG)

    reasonNames = ["impedance" "pop" "bridged" "PREP" "extreme" "muscle" "unknown"];
    headerFields = ["ImpedanceBasedChannelToReject" "PopBasedChannelToReject" "BridgedChannelToReject" ...
        "PREPBasedChannelToReject" "ExtremeDataBasedChannelToReject" "MuscleBasedElectrodesToReject"];

    allLabels = string({EEG.allchan.labels});
    labels = allLabels(~ismember(allLabels, flattenedLabels(EEG.RELAX.ListOfChannelsAfterRejections)));

    baseReasons = repmat("unknown", size(labels));
    for fieldIndex = numel(headerFields):-1:1
        listed = flattenedLabels(rejectionList(EEG, headerFields(fieldIndex)));
        baseReasons(ismember(labels, listed)) = reasonNames(fieldIndex);
    end

    fullReasons = baseReasons;
    [popSubReasonOf, popLabels] = popSubReasons(EEG);
    for labelIndex = find(baseReasons == "pop")
        isPopLabel = popLabels == labels(labelIndex);
        if any(isPopLabel)
            fullReasons(labelIndex) = "pop (" + popSubReasonOf(find(isPopLabel, 1)) + ")";
        end
    end

    reasons = struct('labels', labels, 'reasons', fullReasons, 'summary', "", 'counts', struct(), ...
        'groups', struct('reason', {}, 'labels', {}));
    for reasonName = reasonNames
        reasons.counts.(reasonName) = nnz(baseReasons == reasonName);
        groupReasons = unique(fullReasons(baseReasons == reasonName), 'stable');
        for groupReason = groupReasons
            reasons.groups(end + 1) = struct('reason', groupReason, 'labels', labels(fullReasons == groupReason));
        end
    end
    groupText = strings(1, 0);
    for group = reasons.groups
        groupText(end + 1) = group.reason + ": " + strjoin(group.labels, " ");
    end
    reasons.summary = strjoin(groupText, "; ");
end

function list = rejectionList(EEG, fieldName)
    list = {};
    if isfield(EEG, 'RELAXProcessingExtremeRejections') && isfield(EEG.RELAXProcessingExtremeRejections, fieldName)
        list = EEG.RELAXProcessingExtremeRejections.(fieldName);
    end
end

function [subReasons, popLabels] = popSubReasons(EEG)
    subReasons = strings(1, 0);
    popLabels = strings(1, 0);
    if isfield(EEG.RELAX, 'PopScreen') && isfield(EEG.RELAX.PopScreen, 'badChannels') && isfield(EEG.RELAX.PopScreen, 'badChannelReason')
        popLabels = flattenedLabels(EEG.RELAX.PopScreen.badChannels);
        subReasons = flattenedLabels(EEG.RELAX.PopScreen.badChannelReason);
    end
end

function list = flattenedLabels(value)
    if iscell(value)
        parts = cellfun(@flattenedLabels, value(:)', 'UniformOutput', false);
        list = [strings(1, 0), parts{:}];
    else
        list = string(value);
        list = list(:)';
    end
    list = list(~ismissing(list) & list ~= "");
end
