% RELAX_QC_try_example_windows() - RELAX_QC_example_windows that turns an error into three
%                                 unavailable windows, so one bad example cannot lose a
%                                 long extract run.
%
% Usage:
%   >> windows = RELAX_QC_try_example_windows(pick, summary, paths, cfg, opts);

function windows = RELAX_QC_try_example_windows(pick, summary, paths, cfg, opts)
    try
        windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts);
    catch err
        fprintf('RELAX_QC_try_example_windows: %s failed: %s\n', pick.file, err.message);
        windows = [RELAX_QC_empty_window("blink"), RELAX_QC_empty_window("muscle"), RELAX_QC_empty_window("pop")];
        [windows.note] = deal("Example extraction failed: " + string(err.message));
    end
end
