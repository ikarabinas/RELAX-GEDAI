function A = GEDAI_alpha_run(EEG, refCOV, opts)
% GEDAI_ALPHA_RUN  One broadband-only GEDAI call returning the off, on and removed arms.
%
%   A = GEDAI_alpha_run(EEG, refCOV, opts)
%
%   off = clean + artifacts (the input after GEDAI's reference and high-pass), on = clean,
%   removed = off - on. Also returns threshold and enova (the broadband entries) and the
%   runtime in seconds. Errors if the boundary count or the data size changes.

    p = inputParser;
    addRequired(p, 'EEG', @isstruct);
    addRequired(p, 'refCOV', @isnumeric);
    addRequired(p, 'opts', @isstruct);
    parse(p, EEG, refCOV, opts);

    args = opts.gedaiArgs;
    args{4} = refCOV;
    tic;
    [cleaned, artifacts] = GEDAI(EEG, args{:});
    A.seconds = toc;

    assert(countBoundaries(cleaned) == countBoundaries(EEG), 'GEDAI_alpha:boundaries', ...
        'GEDAI changed the boundary count from %d to %d', countBoundaries(EEG), countBoundaries(cleaned));
    A.on = double(cleaned.data);
    A.removed = double(artifacts.data);
    A.off = A.on + A.removed;
    assert(isequal(size(A.off), size(EEG.data)), 'GEDAI_alpha:size', ...
        'GEDAI returned %s samples for %s input', mat2str(size(A.off)), mat2str(size(EEG.data)));
    A.threshold = cleaned.etc.GEDAI.artifact_threshold_per_band(1);
    A.enova = cleaned.etc.GEDAI.ENOVA_per_band(1);
end

function n = countBoundaries(EEG)
    n = 0;
    if ~isempty(EEG.event)
        n = sum(strcmp({EEG.event.type}, 'boundary'));
    end
end
