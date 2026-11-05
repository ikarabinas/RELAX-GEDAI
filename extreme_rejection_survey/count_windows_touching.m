function n = count_windows_touching(sampleFlags, ranges)
% How many of the [start end] sample ranges contain at least one flagged sample.
    n = 0;
    for i = 1:size(ranges, 1)
        n = n + any(sampleFlags(ranges(i,1):ranges(i,2)));
    end
end
