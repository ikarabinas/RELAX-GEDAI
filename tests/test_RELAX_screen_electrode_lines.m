function test_RELAX_screen_electrode_lines()
% One wrapped list per removal reason, in reason order, then the hole electrodes.
    reasons.groups = struct('reason', {"impedance", "pop (rate)", "muscle"}, 'labels', {["E5" "E6"], "E9", ["E33" "E34" "E35"]});
    lines = RELAX_screen_electrode_lines(reasons, strings(1, 0), 1.4);
    assert(isequal(lines, ["Interpolated for impedance (2): E5, E6"; "Interpolated for pop (rate) (1): E9"; "Interpolated for muscle (3): E33, E34, E35"]));

    withHoles = RELAX_screen_electrode_lines(reasons, ["E5" "E6"], 1.4);
    assert(numel(withHoles) == 4 && withHoles(4) == "Hole electrodes (no real electrode within 1.4× spacing): E5, E6");

    nothing.groups = struct('reason', {}, 'labels', {});
    assert(isempty(RELAX_screen_electrode_lines(nothing, strings(1, 0), 1.4)));

    many.groups = struct('reason', "unknown", 'labels', string(arrayfun(@(k) sprintf('E%d', k), 1:60, 'UniformOutput', false)));
    wrapped = RELAX_screen_electrode_lines(many, strings(1, 0), 1.4);
    assert(numel(wrapped) > 1 && all(strlength(wrapped) <= 110) && startsWith(wrapped(1), "Interpolated for unknown (60):"));
    fprintf('test_RELAX_screen_electrode_lines: PASS\n');
end
