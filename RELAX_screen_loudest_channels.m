% RELAX_screen_loudest_channels() - The channels with the most power inside a span of samples.
%
% Usage:
%   >> labels = RELAX_screen_loudest_channels(filteredData, chanLabels, span, excludeLabels, n);
%
% Inputs:
%   filteredData  - [nChannels x pnts] data, normally the 20-80 Hz band-passed signal
%   chanLabels    - cell of the nChannels labels
%   span          - sample indices (columns of filteredData) to measure over
%   excludeLabels - labels that may not be returned
%   n             - most labels to return
%
% Outputs:
%   labels - [1 x <=n] cell of labels, loudest first, ranked by mean squared value over span

function labels = RELAX_screen_loudest_channels(filteredData, chanLabels, span, excludeLabels, n)
    power = mean(double(filteredData(:, span)).^2, 2);
    candidates = find(~ismember(chanLabels(:), excludeLabels));
    [~, order] = sort(power(candidates), 'descend');
    loudest = candidates(order(1:min(n, numel(order))));
    labels = chanLabels(loudest);
    labels = labels(:)';
end
