function noise = make_pink_noise(nChan, nSamples)
% nChan rows of unit-variance noise whose power falls as 1/f.
    bin = 0:nSamples - 1;
    cyclesPerRecording = max(1, min(bin, nSamples - bin));
    noise = real(ifft(fft(randn(nChan, nSamples), [], 2) ./ sqrt(cyclesPerRecording), [], 2));
    noise = noise ./ std(noise, 0, 2);
end
