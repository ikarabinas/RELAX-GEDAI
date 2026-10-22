% RELAX_screen_reason() - One line saying why a file passed, failed or was not screened.
%
% Usage:
%   >> reason = RELAX_screen_reason(row, meta);
%   >> reason = RELAX_screen_reason(row, meta, baselineNote);
%
% Inputs:
%   row          - one row of the table from RELAX_screen_score_table
%   meta         - its second output (.cutoff, .metricNames)
%   baselineNote - (optional) why the before-wICA trace is missing; appended after " | "
%                  when not empty
%
% Outputs:
%   reason - string. Failing files list each flagging metric with its z; unscreened files
%            list the missing metrics that leave a file unscreened (a missing
%            MuscleStrength never does); passing files state the cutoff.

function reason = RELAX_screen_reason(row, meta, baselineNote)

    if nargin < 3, baselineNote = ""; end
    switch row.ScreenOutcome
        case "notCleaned", reason = "Not cleaned by the pipeline";
        case "fail",       reason = failReason(row, meta);
        case "unscreened", reason = unscreenedReason(row);
        otherwise,         reason = sprintf("All metrics within cutoff %.1f", meta.cutoff);
    end
    if strlength(baselineNote) > 0, reason = reason + " | " + baselineNote; end
end

function reason = failReason(row, meta)
    sentences = strings(1, 0);
    for name = meta.metricNames
        if row.(name + "_flag")
            sentences(end + 1) = sprintf("%s is %.1f robust SD above the dataset median (cutoff %.1f)", ...
                                         name, row.(name + "_z"), meta.cutoff);
        end
    end
    reason = strjoin(sentences, "; ");
end

function reason = unscreenedReason(row)
    metricsThatUnscreen = ["MuscleProportion", "BlinkRatio", "NInterpolated"];
    missing = metricsThatUnscreen(arrayfun(@(name) isnan(row.(name)), metricsThatUnscreen));
    reason = "Not screened";
    if ~isempty(missing), reason = reason + ": " + strjoin(missing, ", ") + " missing"; end
end
