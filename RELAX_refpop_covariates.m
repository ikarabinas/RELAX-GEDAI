% RELAX_refpop_covariates() - Per-file covariates of the reference-pop sub-stage, as
%                             scalars for RELAXProcessingExtremeRejections.
%
% Usage:
%   >> cov = RELAX_refpop_covariates(refInfo, repairInfo);
%
% Every field is a numeric scalar (NaN where nothing was measured) so struct2table can
% append the row. See CLAUDE.md, "Saving covariates for downstream ML classification tasks".

function cov = RELAX_refpop_covariates(refInfo, repairInfo)
    cov = struct('RefPopNCandidates', 0, 'RefPopNSteps', 0, 'RefPopNRejectedUniformity', 0, ...
                 'RefPopNMerged', 0, 'RefPopStepsPerMinute', 0, 'RefPopMedianAmplitude', NaN, 'RefPopMedianRiseMs', NaN, ...
                 'RefPopMedianSignAgreement', NaN, 'RefPopMedianAmplitudeSpread', NaN, ...
                 'RefPopTransitionSeconds', 0, 'RefPopNChannelsUsed', 0, 'RefPopChannelRobustSd', NaN);
    if ~isfield(refInfo, 'available') || ~refInfo.available; return; end
    cov.RefPopNCandidates         = refInfo.nCandidates;
    cov.RefPopNSteps              = refInfo.nSteps;
    cov.RefPopNRejectedUniformity = refInfo.nRejectedUniformity;
    cov.RefPopNMerged             = refInfo.nMerged;
    cov.RefPopStepsPerMinute      = refInfo.stepsPerMinute;
    cov.RefPopNChannelsUsed       = refInfo.channelInfo.nChannelsUsed;
    cov.RefPopChannelRobustSd     = refInfo.channelInfo.robustSd;
    c = refInfo.candidates; acc = logical(c.accepted);
    if any(acc)
        cov.RefPopMedianAmplitude       = median(abs(c.amplitude(acc)));
        cov.RefPopMedianRiseMs          = median(c.riseMs(acc));
        cov.RefPopMedianSignAgreement   = median(c.signAgreement(acc));
        cov.RefPopMedianAmplitudeSpread = median(c.amplitudeSpread(acc));
    end
    if isfield(repairInfo, 'available') && repairInfo.available
        cov.RefPopTransitionSeconds = repairInfo.transitionSeconds;
    end
end