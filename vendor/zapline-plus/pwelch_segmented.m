function [pxx, f] = pwelch_segmented(data, window, noverlap, nfft, fs, boundaries)
% PWELCH_SEGMENTED - Welch PSD that never lets an FFT window straddle a discontinuity.
%
% Averages pwelch estimates computed independently within each stretch of continuous
% data delimited by BOUNDARIES, so the broadband splatter of a discontinuity never
% enters the spectrum.
%
%   data       - samples x channels
%   window     - window vector, or [] to derive one from the segment layout
%   noverlap   - passed through to pwelch ([] = pwelch's default 50%)
%   nfft       - FFT length, scalar or []. MUST NOT be a frequency vector: pwelch
%                treats a vector here as a request for its Goertzel algorithm, which
%                omits the one-sided power doubling and so returns half the PSD.
%                nfft sets the frequency GRID (spacing fs/nfft), which is independent
%                of the WINDOW length that sets the true resolution bandwidth. Callers
%                searching a narrow frequency window must pass an nfft fine enough
%                that the window actually contains bins.
%   fs         - sampling rate in Hz
%   boundaries - sample indices of discontinuities, or [] for continuous data
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    if isempty(window)
        if isempty(boundaries)
            window = hanning(size(data, 1));
        else
            bounds_vec = boundaries(:)';
            segment_lengths = diff(unique(sort([1, bounds_vec, size(data, 1) + 1])));
            max_win_length = max(1, min(segment_lengths));
            window = hanning(max_win_length);
        end
    end

    % pwelch wraps the data modulo nfft when nfft is shorter than the window, which
    % silently coarsens the grid instead of zero-padding. Never allow that.
    if isempty(nfft)
        nfft = 2^nextpow2(length(window));
    else
        nfft = max(nfft, 2^nextpow2(length(window)));
    end

    if isempty(boundaries) || length(boundaries) < 2
        [pxx, f] = pwelch(data, window, noverlap, nfft, fs);
        return;
    end

    boundaries = unique(sort([1; boundaries(:); (size(data, 1) + 1)]));
    pxx_total = 0;
    total_segments = 0;
    f = [];

    for i = 1:length(boundaries) - 1
        start_idx = boundaries(i);
        end_idx   = boundaries(i + 1) - 1;
        seg_len   = end_idx - start_idx + 1;

        if seg_len >= length(window)
            % NOTE: nfft is passed in unchanged every iteration, and the frequency
            % vector comes back into a separate variable. Assigning f into the nfft
            % argument would switch pwelch to Goertzel from the second segment on.
            [pxx_seg, f] = pwelch(data(start_idx:end_idx, :), window, noverlap, nfft, fs);
            pxx_total = pxx_total + pxx_seg;
            total_segments = total_segments + 1;
        end
    end

    if total_segments > 0
        pxx = pxx_total / total_segments;
    else
        [pxx, f] = pwelch(data, window, noverlap, nfft, fs);
    end
end
