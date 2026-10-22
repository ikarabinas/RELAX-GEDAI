% RELAX_screen_outcome() - Pass/fail/unscreened/notCleaned verdict per file from its metric z-scores.
%
% Usage:
%   >> [outcome, flaggedBy] = RELAX_screen_outcome(pipelineOutcome, values, z, cutoff, ...
%                                                  metricNames, missingMakesUnscreened);
%
% Inputs:
%   pipelineOutcome        - [n x 1] string, "cleaned" for files the pipeline cleaned
%   values, z              - [n x m] metric values and their robust z-scores
%   cutoff                 - a metric flags when z > cutoff (strictly)
%   metricNames            - [1 x m] string
%   missingMakesUnscreened - [1 x m] logical, whether a NaN value in that metric leaves
%                            a file that no metric flags unscreened
%
% Outputs:
%   outcome   - [n x 1] string; precedence "notCleaned", "fail", "unscreened", "pass"
%   flaggedBy - [n x 1] string, flagging metric names joined with "; " ("" when none)

function [outcome, flaggedBy] = RELAX_screen_outcome(pipelineOutcome, values, z, cutoff, metricNames, missingMakesUnscreened)

    nFiles = numel(pipelineOutcome);
    notCleaned = pipelineOutcome ~= "cleaned";
    flagged = z > cutoff & ~notCleaned;
    failed = any(flagged, 2);
    unscreened = any(isnan(values) & missingMakesUnscreened, 2);

    outcome = repmat("pass", nFiles, 1);
    outcome(unscreened) = "unscreened";
    outcome(failed) = "fail";
    outcome(notCleaned) = "notCleaned";

    flaggedBy = strings(nFiles, 1);
    for f = find(failed)'
        flaggedBy(f) = strjoin(metricNames(flagged(f, :)), "; ");
    end
end
