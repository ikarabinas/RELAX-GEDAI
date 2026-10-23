% RELAX_screen_muscle_caption() - The caption under one muscle window of a review sheet.
%
% Usage:
%   >> caption = RELAX_screen_muscle_caption(k, startSeconds, slope, channel, threshold);
%
% Inputs:
%   k            - the window's number on the sheet
%   startSeconds - start of the window's epoch in the continuous data
%   slope        - the highest muscle slope over the epoch's channels
%   channel      - the channel holding that slope
%   threshold    - cfg.MuscleSlopeThreshold
%
% A slope that is not above the threshold is captioned ", below threshold".

function caption = RELAX_screen_muscle_caption(k, startSeconds, slope, channel, threshold)

    caption = sprintf("Muscle epoch %d, starts at %.1f s: slope %.2f on %s (threshold %.2f)", k, startSeconds, slope, channel, threshold);
    if ~(slope > threshold)
        caption = caption + ", below threshold";
    end
end
