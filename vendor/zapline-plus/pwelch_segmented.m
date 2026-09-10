function [pxx, f] = pwelch_segmented(data, window, noverlap, f, fs, boundaries)
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

    if isempty(boundaries) || length(boundaries) < 2
        [pxx, f] = pwelch(data, window, noverlap, f, fs);
        return;
    end
    
    boundaries = unique(sort([1; boundaries(:); (size(data, 1)+1)]));
    pxx_total = 0;
    total_chunks = 0;

    for i = 1:length(boundaries)-1
        start_idx = boundaries(i);
        end_idx = boundaries(i+1) - 1;
        chunk_len = end_idx - start_idx + 1;
        
        if chunk_len >= length(window)
            [pxx_chunk, f] = pwelch(data(start_idx:end_idx, :), window, noverlap, f, fs);
            pxx_total = pxx_total + pxx_chunk;
            total_chunks = total_chunks + 1;
        end
    end
    
    if total_chunks > 0
        pxx = pxx_total / total_chunks;
    else
        [pxx, f] = pwelch(data, window, noverlap, f, fs);
    end
end
