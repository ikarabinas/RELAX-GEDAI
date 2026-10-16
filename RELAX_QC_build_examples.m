% RELAX_QC_build_examples() - Figure 1 example windows for every picked file.
%
% Usage:
%   >> examples = RELAX_QC_build_examples(picks, T, summaries, cfgs, processed, inputs, opts);
%
% T, summaries and cfgs are the cache tables, its summaries (row-matched to T.File) and the
% participant configs from RELAX_QC_load_tables. A pick whose summary is not "ok" is
% warned about and left out.

function examples = RELAX_QC_build_examples(picks, T, summaries, cfgs, processed, inputs, opts)
    examples = struct('file', {}, 'role', {}, 'windows', {});
    for k = 1:numel(picks)
        i = find(T.File == picks(k).file, 1);
        if summaries(i).status ~= "ok"
            warning('RELAX_QC:exampleUnavailable', 'Example %s has status %s; not shown.', picks(k).file, summaries(i).status);
            continue
        end
        windows = RELAX_QC_try_example_windows(picks(k), summaries(i), RELAX_QC_file_paths(processed, inputs, T.File(i)), ...
                                               cfgs(char(T.Participant(i))), opts);
        examples(end+1) = struct('file', picks(k).file, 'role', picks(k).role, 'windows', windows); %#ok<AGROW>
        reportFidelity(picks(k).file, windows);
    end
end

function reportFidelity(file, windows)
    for w = windows([windows.kind] == "pop" | [windows.kind] == "popKept")
        if ~w.available, fprintf('%s %s: %s\n', file, w.kind, w.note); continue, end
        fprintf('%s %s repaired %s: relative RMS vs Extremes_Rejected %.3f (all shown channels: %s); %.0f%% of window compared\n', ...
                file, w.kind, w.channels{1}, w.fidelity.relativeRms(1), mat2str(w.fidelity.relativeRms, 3), 100 * w.fidelity.fractionCompared);
    end
end
