function test_RELAX_screen_tick_labels()
% Only interpolated channels are marked, each with its own reason.
    reasons = struct('labels', ["E5" "E37"], 'reasons', ["impedance" "pop (rate)"]);
    labels = RELAX_screen_tick_labels({'E1'; 'E37'; 'E5'}, reasons);
    assert(isequal(labels, {'E1'; 'E37 (interp., pop (rate))'; 'E5 (interp., impedance)'}));
    none = struct('labels', strings(1, 0), 'reasons', strings(1, 0));
    assert(isequal(RELAX_screen_tick_labels({'E1'; 'E2'}, none), {'E1'; 'E2'}));
    fprintf('test_RELAX_screen_tick_labels: PASS\n');
end
