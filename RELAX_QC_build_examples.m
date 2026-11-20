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
        reportUnavailablePops(picks(k).file, windows);
    end
end

function reportUnavailablePops(file, windows)
    for w = windows([windows.kind] == "pop" | [windows.kind] == "popKept")
        if ~w.available, fprintf('%s %s: %s\n', file, w.kind, w.note); end
    end
end
