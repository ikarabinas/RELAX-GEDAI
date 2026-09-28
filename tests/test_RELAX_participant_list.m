function test_RELAX_participant_list()
% RELAX_SET_PARAMETERS_AND_RUN_slurm.m maps SLURM array index -> participant through a
% cached list (participant_list.mat) so concurrent array tasks agree on the mapping. The
% cache used to be written once and trusted forever: when m579's files were added to a
% folder whose cache held 25 IDs, array task 26 printed "no participant to process" and
% m579 was silently never cleaned (batch 19739567). The helper must keep every cached
% index where it is and append newcomers at the end, so an index never changes meaning
% between runs while a grown folder still gets its new participants.

    tmp = tempname(); mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's'));
    cacheFile = fullfile(tmp, 'participant_list.mat');

    % 1. no cache: stable order of first appearance, cache written
    ids = {'m20','m20','m05','m11','m05'};
    lastwarn('');
    list = RELAX_participant_list(cacheFile, ids);
    assert(isequal(list, {'m20','m05','m11'}), 'stable first-appearance order');
    assert(isrow(list), 'row cell');
    assert(exist(cacheFile, 'file') == 2, 'cache written');
    assert(isempty(lastwarn()), 'no warning on first creation');

    % 2. unchanged folder: same list, no warning
    lastwarn('');
    list = RELAX_participant_list(cacheFile, ids);
    assert(isequal(list, {'m20','m05','m11'}));
    assert(isempty(lastwarn()), 'no warning when nothing changed');

    % 3. a new participant that sorts before the others: appended, indices preserved
    lastwarn('');
    list = RELAX_participant_list(cacheFile, [{'m01'} ids {'m01'}]);
    assert(isequal(list, {'m20','m05','m11','m01'}), 'newcomer appended at the end, got %s', strjoin(list, ','));
    [~, wid] = lastwarn();
    assert(strcmp(wid, 'RELAX_participant_list:added'), 'warning names the addition');
    s = load(cacheFile, 'uniqueParticipants');
    assert(isequal(s.uniqueParticipants(:)', {'m20','m05','m11','m01'}), 'cache updated with the appended list');

    % 4. a cached participant whose files are gone: index kept (so the others do not
    %    shift), warning issued, cache untouched
    lastwarn('');
    list = RELAX_participant_list(cacheFile, {'m20','m11','m01'});
    assert(isequal(list, {'m20','m05','m11','m01'}), 'missing participant keeps its slot');
    [~, wid] = lastwarn();
    assert(strcmp(wid, 'RELAX_participant_list:missing'), 'warning names the missing participant');

    % 5. column input works the same
    list = RELAX_participant_list(cacheFile, {'m20';'m05';'m11';'m01'});
    assert(isequal(list, {'m20','m05','m11','m01'}));

    fprintf('test_RELAX_participant_list: PASS\n');
end
