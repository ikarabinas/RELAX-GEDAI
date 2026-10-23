function test_RELAX_screen_channel_reasons()
    r = RELAX_screen_channel_reasons(header());
    assert(isequal(r.labels, ["E2" "E3" "E5" "E6" "E7"]));
    assert(isequal(r.reasons, ["pop (rate)" "pop (flat)" "impedance" "muscle" "unknown"]));
    assert(r.summary == "impedance: E5; pop (rate): E2; pop (flat): E3; muscle: E6; unknown: E7");
    assert(r.counts.pop == 2 && r.counts.unknown == 1 && r.counts.bridged == 0);
    assert(r.counts.impedance == 1 && r.counts.PREP == 0 && r.counts.extreme == 0 && r.counts.muscle == 1);
    assert(numel(r.groups) == 5);
    assert(isequal([r.groups.reason], ["impedance" "pop (rate)" "pop (flat)" "muscle" "unknown"]));
    assert(isequal(r.groups(2).labels, "E2") && isequal(r.groups(4).labels, "E6"));

    noRejections = rmfield(header(), 'RELAXProcessingExtremeRejections');
    r = RELAX_screen_channel_reasons(noRejections);
    assert(isequal(r.reasons, repmat("unknown", 1, 5)) && r.counts.unknown == 5 && r.counts.pop == 0);

    nothing = header();
    nothing.RELAX.ListOfChannelsAfterRejections = {nothing.allchan.labels};
    r = RELAX_screen_channel_reasons(nothing);
    assert(isempty(r.labels) && isempty(r.reasons) && r.summary == "" && isempty(r.groups));
    assert(all(structfun(@(count) count == 0, r.counts)));

    noPopScreen = header();
    noPopScreen.RELAX = rmfield(noPopScreen.RELAX, 'PopScreen');
    r = RELAX_screen_channel_reasons(noPopScreen);
    assert(isequal(r.reasons, ["pop" "pop" "impedance" "muscle" "unknown"]));
    assert(r.summary == "impedance: E5; pop: E2 E3; muscle: E6; unknown: E7");

    flatFirst = header();
    flatFirst.RELAX.PopScreen.badChannels = {'E2'; 'E3'};
    flatFirst.RELAX.PopScreen.badChannelReason = {'flat'; 'rate'};
    r = RELAX_screen_channel_reasons(flatFirst);
    assert(isequal([r.groups.reason], ["impedance" "pop (flat)" "pop (rate)" "muscle" "unknown"]), 'pop sub-reasons follow the first label in allchan');

    bareFields = header();
    bareFields.RELAXProcessingExtremeRejections = struct('MuscleBasedElectrodesToReject', 'E6', 'PREPBasedChannelToReject', "E7");
    r = RELAX_screen_channel_reasons(bareFields);
    assert(isequal(r.reasons, ["unknown" "unknown" "unknown" "muscle" "PREP"]), 'absent fields are empty lists; char and string lists are read');

    twoScreens = header();
    twoScreens.RELAXProcessingExtremeRejections.ExtremeDataBasedChannelToReject = {'E5'};
    r = RELAX_screen_channel_reasons(twoScreens);
    assert(r.reasons(3) == "impedance" && r.counts.extreme == 0, 'a channel in two lists takes the earlier reason');
    fprintf('test_RELAX_screen_channel_reasons: PASS\n');
end

function EEG = header()
    labels = compose('E%d', (1:8)');
    EEG.allchan = struct('labels', labels);
    EEG.RELAX.ListOfChannelsAfterRejections = {'E1'; 'E4'; 'E8'};
    EEG.RELAX.PopScreen.badChannels = {'E3'; 'E2'};
    EEG.RELAX.PopScreen.badChannelReason = {'flat'; 'rate'};
    EEG.RELAXProcessingExtremeRejections.ImpedanceBasedChannelToReject = {{'E5'}};
    EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject = {{'E2'; 'E3'}};
    EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject = {{}};
    EEG.RELAXProcessingExtremeRejections.ExtremeDataBasedChannelToReject = {};
    EEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject = {'E6'};
end
