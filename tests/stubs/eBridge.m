function EB = eBridge(EPO, varargin)
% Test double for the EEGLAB ebridge plugin's eBridge(), used only by
% tests/test_RELAX_bridged_channels.m so that RELAX_bridged_channels can be
% exercised end-to-end without the real ebridge plugin installed. Branches on the
% channel count of EPO.data to return a canned result per test scenario - see
% tests/test_RELAX_bridged_channels.m for what each channel count means.

    nChan = size(EPO.data, 1);
    labels = {EPO.chanlocs.labels};

    switch nChan
        case 40 % "no bridging": no near-zero peak in the ED distribution
            EB.Info.EDcutoff = -1;
            EB.Info.EDscale = NaN;
            EB.Bridged.Count = 0;
            EB.Bridged.Pairs = cell(2,0);
        case 41 % "clean file": a cutoff was found but nothing crossed it
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.0;
            EB.Bridged.Count = 0;
            EB.Bridged.Pairs = cell(2,0);
        case 30 % a single 2-channel bridge, well under any QC cap
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.2;
            EB.Bridged.Pairs = {labels{1}; labels{2}};
            EB.Bridged.Count = 1;
        case 12 % a 5-channel chain: triggers both oversizedCluster and excessiveBridging
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.2;
            EB.Bridged.Pairs = [labels(1:4); labels(2:5)];
            EB.Bridged.Count = 4;
        otherwise
            error('eBridge test double: unexpected channel count %d - add a case for it', nChan);
    end
end
