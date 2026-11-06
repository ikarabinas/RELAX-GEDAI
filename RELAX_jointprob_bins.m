function b = RELAX_jointprob_bins(x, lo, hi, bins)
% The histogram bin EEGLAB's realproba assigns each value of x for a channel spanning
% [lo, hi]. Values outside that span fall outside 1..bins.
    b = floor((x - lo) / (hi - lo) * (bins - 1)) + 1;
end
