% RELAX_QC_fig_paired_metrics() - Figure 3: blink and muscle metrics before and after
%                                 wICA, one line per file.

function caption = RELAX_QC_fig_paired_metrics(cache, outDir)

    style = RELAX_QC_style();
    ok = cache.summaries([cache.summaries.status] == "ok");
    panels = {'blinkRatio', 'Blink amplitude ratio', 1; ...
              'muscleProportion', {'Epochs with muscle', 'above threshold (%)'}, 100; ...
              'muscleStrength', {'Muscle slope', 'above threshold'}, 1};

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 18 9], 'Color', 'w');
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:size(panels, 1)
        ax = nexttile(layout);
        if isempty(ok)
            RELAX_QC_draw_note(ax, "No aligned files");
            continue
        end
        hold(ax, 'on');
        before = arrayfun(@(s) s.metricsBaseline.(panels{p, 1}), ok) * panels{p, 3};
        after  = arrayfun(@(s) s.metricsCleaned.(panels{p, 1}), ok) * panels{p, 3};
        plot(ax, [1 2], [before(:) after(:)]', 'Color', [style.before 0.5], 'LineWidth', style.thin);
        plot(ax, [1 2], [median(before, 'omitnan') median(after, 'omitnan')], '-o', 'Color', style.cleaned, ...
             'LineWidth', style.bold, 'MarkerSize', 6, 'MarkerFaceColor', style.cleaned);
        xlim(ax, [0.7 2.3]); xticks(ax, [1 2]); xticklabels(ax, {'Before wICA', 'Cleaned'}); xtickangle(ax, 0);
        ylabel(ax, panels{p, 2}); ax.YGrid = 'on'; ax.FontSize = style.fontSize;
    end
    RELAX_QC_save_figure(fig, outDir, 'fig3_paired_metrics');

    caption = string(sprintf(['## Figure 3 — Artifact metrics before and after wICA\n\n' ...
        'Each grey line is one file (%d files); the blue line joins the medians. "Before wICA" is the ' ...
        'sample-matched baseline, which is after pop repair, ZaplinePlus and extreme-period excision, not raw data. ' ...
        'Both sides are computed by the pipeline''s own metric functions. Blink amplitude ratio: mean absolute ' ...
        'amplitude within ±0.5 s of each blink peak over that in the 0.5 s at either end of a ±2 s epoch, ' ...
        'averaged over the blink electrodes; 1 means no blink is left. Muscle: share of epochs whose ' ...
        'log-frequency/log-power slope exceeds MuscleSlopeThreshold (%.2f) in any channel, and the mean amount by ' ...
        'which super-threshold slopes exceed it. Good cleaning: every line falls, blink ratio ends near 1.\n'], ...
        numel(ok), cache.meta.muscleSlopeThreshold));
end
