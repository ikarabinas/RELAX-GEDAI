function uniqueParticipants = RELAX_participant_list(cacheFile, participantIDs)
% RELAX_participant_list() - SLURM array index -> participant mapping, pinned by a cache.
%
% Usage:
%   >> uniqueParticipants = RELAX_participant_list(cacheFile, participantIDs);
%
% Inputs:
%   cacheFile      - full path of the .mat file that pins the mapping (variable
%                    'uniqueParticipants'); created if absent
%   participantIDs - cell array of the participant ID parsed from every input file
%                    (one entry per file, duplicates expected)
%
% Output:
%   uniqueParticipants - 1 x N cell row; array task k processes uniqueParticipants{k}
%
% Every concurrently running array task must agree on which index is which
% participant, and a later rerun of "--array=7" must mean the same participant it
% meant last time, so the mapping is cached rather than recomputed. The cache used to
% be written once and trusted forever, which silently dropped any participant whose
% files arrived after the first run (m579, batch 19739567: task 26 found "no
% participant to process"). Now:
%   - a participant already in the cache keeps its index, always;
%   - a participant present in the folder but not in the cache is appended at the end
%     and the cache is rewritten (warning RELAX_participant_list:added - the sbatch
%     --array range must cover the new length);
%   - a cached participant with no files left keeps its slot so nobody else shifts;
%     its task will find nothing to do (warning RELAX_participant_list:missing).
% The rewrite goes through a temp file + movefile so a task starting mid-write never
% loads a truncated cache.

    current = unique(participantIDs(:)', 'stable');

    if exist(cacheFile, 'file') == 2
        s = load(cacheFile, 'uniqueParticipants');
        cached = s.uniqueParticipants(:)';
        added = current(~ismember(current, cached));
        missing = cached(~ismember(cached, current));
        uniqueParticipants = [cached added];
        if ~isempty(added)
            warning('RELAX_participant_list:added', ...
                ['%d participant(s) not in %s appended at index %d-%d: %s. Make sure the ' ...
                 'sbatch --array range reaches %d.'], numel(added), cacheFile, ...
                numel(cached) + 1, numel(uniqueParticipants), strjoin(added, ', '), ...
                numel(uniqueParticipants));
            local_save(cacheFile, uniqueParticipants);
        end
        if ~isempty(missing)
            warning('RELAX_participant_list:missing', ...
                ['%d cached participant(s) have no files in the folder any more: %s. Their ' ...
                 'index is kept so the others do not shift; the corresponding array task(s) ' ...
                 'will have nothing to process.'], numel(missing), strjoin(missing, ', '));
        end
    else
        uniqueParticipants = current;
        local_save(cacheFile, uniqueParticipants);
    end
end

function local_save(cacheFile, uniqueParticipants)
    tmpFile = [tempname(fileparts(cacheFile)) '.mat'];
    save(tmpFile, 'uniqueParticipants');
    movefile(tmpFile, cacheFile, 'f');
end
