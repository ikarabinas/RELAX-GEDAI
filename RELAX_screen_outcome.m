% RELAX_screen_outcome() - Pass/fail/unscreened/notCleaned verdict per file from its metric flags.
%
% Usage:
%   >> [outcome, flaggedBy] = RELAX_screen_outcome(pipelineOutcome, values, flagged, ...
%                                                  metricNames, missingMakesUnscreened);
%
% Inputs:
%   pipelineOutcome        - [n x 1] string, "cleaned" for files the pipeline cleaned
%   values                 - [n x m] metric values
%   flagged                - [n x m] logical, whether each metric flags that file (rows of files
%                            that were not cleaned are ignored)
%   metricNames            - [1 x m] string
%   missingMakesUnscreened - [1 x m] logical, whether a NaN value in that metric leaves
%                            a file that no metric flags unscreened
%
% Outputs:
%   outcome   - [n x 1] string; precedence "notCleaned", "fail", "unscreened", "pass"
%   flaggedBy - [n x 1] string, flagging metric names joined with "; " ("" when none)

function [outcome, flaggedBy] = RELAX_screen_outcome(pipelineOutcome, values, flagged, metricNames, missingMakesUnscreened)

    nFiles = numel(pipelineOutcome);
    notCleaned = pipelineOutcome ~= "cleaned";
    flagged = flagged & ~notCleaned;
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
