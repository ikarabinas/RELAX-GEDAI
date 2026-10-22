% RELAX_QC_fig_spectra() - Figure 2: dataset power spectra of the input, the pre-wICA
%                          baseline, the cleaned data and what wICA removed.
%
% Usage:
%   >> caption = RELAX_QC_fig_spectra(cache, outDir);

function caption = RELAX_QC_fig_spectra(cache, outDir)

    style = RELAX_QC_style();
    ok = cache.summaries([cache.summaries.status] == "ok");
    caption = string(sprintf(['## Figure 2 — Power spectra across the dataset\n\n' ...
        'Median across %d files (bold) with interquartile band and one faint line per file; power is averaged ' ...
        'over channels, from 2 s Hann windows that never span a boundary. The input spectrum is from the unprocessed ' ...
        'file and is not sample-aligned with the others. Good cleaning: the removed spectrum looks like artifact ' ...
        '(low-frequency eye power, a broadband rise above about 20 Hz) with no alpha peak, and the right panel stays ' ...
        'near 0 dB in the shaded 8–12 Hz band.\n'], numel(ok)));

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 20 9], 'Color', 'w');
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    if isempty(ok)
        RELAX_QC_draw_note(nexttile(layout, [1 2]), "No aligned files");
    else
        freqs = ok(1).freqs;
        band = freqs >= 1 & freqs <= 80;
        f = freqs(band);
        series = {'Input (not sample-aligned)', 'psdInput', style.input; ...
                  'Before wICA', 'psdBaseline', style.before; ...
                  'Cleaned', 'psdCleaned', style.cleaned; ...
                  'Removed by wICA + GEDAI', 'psdRemoved', style.removed};

        ax = nexttile(layout); hold(ax, 'on');
        for s = 1:size(series, 1)
            P = log10(vertcat(ok.(series{s, 2})));
            if all(isnan(P(:, band)), 'all'), continue, end
            medianBand(ax, f, P(:, band), series{s, 3}, series{s, 1});
        end
        xlabel(ax, 'Frequency (Hz)'); ylabel(ax, {'log_{10} power', '(\muV^2/Hz)'});
        legend(ax, 'Location', 'northoutside', 'NumColumns', 2, 'Box', 'off'); title(ax, 'Channel-mean spectra');

        ax = nexttile(layout); hold(ax, 'on');
        ratio = 10 * log10(vertcat(ok.psdCleaned) ./ vertcat(ok.psdBaseline));
        lims = quantile(ratio(:, band), [0.01 0.99], 'all');
        yLimits = [min(-1, lims(1)) max(1, lims(2))];
        patch(ax, [8 12 12 8], [-1e3 -1e3 1e3 1e3], style.muted, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        yline(ax, 0, '-', 'Color', style.muted, 'HandleVisibility', 'off');
        medianBand(ax, f, ratio(:, band), style.cleaned, 'Cleaned / before wICA');
        ylim(ax, yLimits);
        xlabel(ax, 'Frequency (Hz)'); ylabel(ax, {'Power change', '(dB)'}); title(ax, 'Cleaned relative to before wICA');
        legend(ax, 'Location', 'northoutside', 'Box', 'off');

        set(findobj(fig, 'Type', 'axes'), 'XLim', [1 80], 'FontSize', style.fontSize, 'YGrid', 'on');
    end
    RELAX_QC_save_figure(fig, outDir, 'fig2_spectra');
end

function medianBand(ax, f, P, colour, name)
    style = RELAX_QC_style();
    plot(ax, f, P', 'Color', [colour 0.12], 'LineWidth', 0.4, 'HandleVisibility', 'off');
    q = quantile(P, [0.25 0.5 0.75], 1);
    if size(P, 1) == 1, q = [P; P; P]; end
    fill(ax, [f fliplr(f)], [q(1, :) fliplr(q(3, :))], colour, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, f, q(2, :), 'Color', colour, 'LineWidth', style.bold, 'DisplayName', name);
end
