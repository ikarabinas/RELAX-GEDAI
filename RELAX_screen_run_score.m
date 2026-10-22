% RELAX_screen_run_score() - Score stage of the residual-noise screen: measure every file,
%                            write residual_screen.mat/.csv and a README.
%
% Usage:
%   >> [screen, meta] = RELAX_screen_run_score(opts);
%
% opts: built by RELAX_residual_screen (.myPath, .outputFolderName, .outDir, .cutoff,
% .sampleSize, .seed). A cleaned file whose header cannot be read keeps its muscle metrics,
% gets NaN BlinkRatio and NInterpolated (so it is unscreened), and the error message goes in
% SheetMessage. See RELAX_screen_score_table for the columns of screen and meta.

function [screen, meta] = RELAX_screen_run_score(opts)

    processed = fullfile(opts.myPath, opts.outputFolderName);
    [T, cfgs] = RELAX_QC_load_tables(fullfile(processed, 'PerParticipantStats'), fullfile(processed, 'Cleaned_Data'));
    cleanedDir = fullfile(processed, 'Cleaned_Data');

    nFiles = height(T);
    metrics = repmat(nanMetrics(), nFiles, 1);
    headerMessages = strings(nFiles, 1);
    for f = find(T.Outcome == "cleaned")'
        cleanedPath = fullfile(cleanedDir, T.File(f) + RELAX_QC_cleaned_suffix(cleanedDir, T.File(f)));
        [metrics(f), headerMessages(f)] = measureOneFile(cleanedPath, T(f, :), cfgs(char(T.Participant(f))));
        fprintf('RELAX_screen_run_score %d/%d %s %s\n', f, nFiles, T.File(f), headerMessages(f));
    end

    [screen, meta] = RELAX_screen_score_table(T, metrics, opts);
    screen.SheetMessage = headerMessages;

    if ~isfolder(opts.outDir), mkdir(opts.outDir); end
    save(fullfile(opts.outDir, 'residual_screen.mat'), 'screen', 'meta');
    writeCsv(screen, fullfile(opts.outDir, 'residual_screen.csv'));
    writeReadme(screen, meta, fullfile(opts.outDir, 'README.md'));
    fprintf('RELAX_screen_run_score: %s\n', outcomeCounts(screen));
end

function m = nanMetrics()
    m = struct('MuscleProportion', NaN, 'MuscleStrength', NaN, 'BlinkRatio', NaN, 'NInterpolated', NaN);
end

function [m, message] = measureOneFile(cleanedPath, row, cfg)
    [header, message] = loadHeader(cleanedPath);
    if message == ""
        m = RELAX_screen_file_metrics(row, cfg, {header.allchan.labels}, ...
                                      string(header.RELAX.ListOfChannelsAfterRejections));
    else
        m = nanMetrics();
        m.MuscleProportion = row.TableCleanedMuscleProportion;
        m.MuscleStrength = row.TableCleanedMuscleStrength;
        warning('RELAX_screen:headerUnreadable', '%s: %s', row.File, message);
    end
end

function [header, message] = loadHeader(cleanedPath)
    header = [];
    message = "";
    try
        header = pop_loadset('filename', char(cleanedPath), 'loadmode', 'info');
    catch err
        message = string(err.message);
    end
end

function writeCsv(screen, csvPath)
    forCsv = screen;
    for column = string(forCsv.Properties.VariableNames)
        if islogical(forCsv.(column)), forCsv.(column) = double(forCsv.(column)); end
    end
    writetable(forCsv, csvPath);
end

function writeReadme(screen, meta, readmePath)
    lines = [ ...
        "# Residual-noise screen", ...
        "", ...
        sprintf("Run on %s with cutoff %g, seed %d.", string(meta.created, 'yyyy-MM-dd'), meta.cutoff, meta.seed), ...
        sprintf("%d cleaned files were scored; the median and scale of each metric come from those files only.", meta.nCleaned), ...
        "", ...
        "## Outcomes", ...
        "", ...
        "- " + outcomeCounts(screen), ...
        "", ...
        "A metric flags a file when its robust z exceeds the cutoff; the <metric>_flag columns hold those flags. " + ...
        "A file is unscreened when a metric other than the muscle slope is missing and no metric flags it.", ...
        "", ...
        "The <metric>_scale columns are the robust SD that z divides by: 1.4826 * MAD, or 1.2533 * mean|x - median| when MAD is 0.", ...
        ""];
    for metric = meta.zeroScaleMetrics
        lines(end + 1) = sprintf("%s has no spread across the cleaned files, so it cannot flag any file in this run.", metric); %#ok<AGROW>
    end
    lines = [lines, "", ...
        "## Reading a sheet", ...
        "", ...
        "Each review sheet stacks 3 muscle windows and 3 blink windows of 3 s, taken from the cleaned data. " + ...
        "Blue traces are cleaned data and grey traces are the data before wICA. " + ...
        "The sheet names the channels shown and the topography of the strongest blink. " + ...
        "A file is acceptable when the windows show no residual muscle bursts or blink deflections that stand out from the surrounding EEG."];
    writelines(lines, readmePath);
end

function text = outcomeCounts(screen)
    outcomes = ["pass", "fail", "unscreened", "notCleaned"];
    counts = arrayfun(@(o) sum(screen.ScreenOutcome == o), outcomes);
    text = strjoin(outcomes + " " + counts, ", ");
end
