function test_RELAX_screen_file_metrics()
    cfg = struct('BlinkElectrodes', {{'E1'; 'E2'; 'E3'}});
    allchan = {'E1', 'E2', 'E3', 'E12', 'E50'};
    row = table(0.2, -0.25, {[1.1; NaN; 1.3; 2.0; 9]}, 0, 'VariableNames', ...
        {'TableCleanedMuscleProportion', 'TableCleanedMuscleStrength', 'TableCleanedBlinkRatio', 'BlinkElectrodeFallback'});
    [m, note] = RELAX_screen_file_metrics(row, cfg, allchan, {'E1', 'E3', 'E50'});
    assert(m.MuscleProportion == 0.2 && m.MuscleStrength == -0.25);
    assert(abs(m.BlinkRatio - 1.2) < 1e-12 && m.NInterpolated == 2 && note == "");

    row.BlinkElectrodeFallback = 1;
    m = RELAX_screen_file_metrics(row, cfg, allchan, {'E1', 'E3', 'E50'});
    assert(abs(m.BlinkRatio - mean([1.1 1.3 2.0])) < 1e-12, 'fallback row adds E12');

    row.BlinkElectrodeFallback = 2;
    [m, note] = RELAX_screen_file_metrics(row, cfg, allchan, allchan);
    assert(isnan(m.BlinkRatio) && contains(note, "no blink electrodes"));

    row.BlinkElectrodeFallback = 0;
    row.TableCleanedBlinkRatio = {[1.1; 1.2]};
    [m, note] = RELAX_screen_file_metrics(row, cfg, allchan, allchan);
    assert(isnan(m.BlinkRatio) && note == "Blink-ratio vector has 2 entries, allchan has 5", 'length mismatch is explained');

    row.BlinkElectrodeFallback = NaN; row.TableCleanedBlinkRatio = {NaN};
    m = RELAX_screen_file_metrics(row, cfg, allchan, allchan);
    assert(isnan(m.BlinkRatio) && m.NInterpolated == 0);
    [~, note] = RELAX_screen_file_metrics(row, cfg, allchan, allchan);
    assert(note == "No cleaned blink ratio in the tables");
    assert(isequal(RELAX_backup_blink_electrodes(), {'E12'; 'E20'; 'E27'; 'E34'; 'E21'}));

    configured = {'E1'; 'E2'; 'E3'};
    assert(isequal(RELAX_screen_blink_labels(cfg, 0), configured));
    assert(isequal(RELAX_screen_blink_labels(cfg, NaN), configured));
    assert(isequal(RELAX_screen_blink_labels(cfg, 1), [configured; RELAX_backup_blink_electrodes()]));
    assert(isempty(RELAX_screen_blink_labels(cfg, 2)));
    fprintf('test_RELAX_screen_file_metrics: PASS\n');
end
