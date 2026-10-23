function test_RELAX_screen_hole_depth()
% Flat triangular lattice with spacing 1: E1 sits at the centre, with its six neighbours at distance 1.
    locs = lattice();
    ring = labelsWithin(locs, 1);

    hole = RELAX_screen_hole_depth(locs, {}, 1.4);
    assert(hole.depth == 0 && hole.deepest == "" && isempty(hole.electrodes) && abs(hole.spacing - 1) < 1e-9 && hole.note == "");

    hole = RELAX_screen_hole_depth(locs, {'E1'}, 1.4);
    assert(abs(hole.depth - 1) < 1e-9 && hole.deepest == "E1" && isempty(hole.electrodes));

    hole = RELAX_screen_hole_depth(locs, rowOfFive(locs), 1.4);
    assert(abs(hole.depth - 1) < 1e-9 && isempty(hole.electrodes), 'a line never leaves an electrode without a real neighbour');

    hole = RELAX_screen_hole_depth(locs, ring, 1.4);
    assert(abs(hole.depth - sqrt(3)) < 1e-9 && hole.deepest == "E1" && isequal(hole.electrodes, "E1"), 'centre plus ring');

    noX = locs; noX(strcmp({noX.labels}, 'E2')).X = [];
    hole = RELAX_screen_hole_depth(noX, {'E1'}, 1.4);
    assert(abs(hole.depth - 1) < 1e-9, 'a real electrode without coordinates is skipped, and nothing after it shifts');

    hole = RELAX_screen_hole_depth(noX, {'E2'}, 1.4);
    assert(isnan(hole.depth) && hole.note == "Interpolated electrodes without coordinates: E2");

    bare = locs(1:3); [bare(2:3).X] = deal([]);
    hole = RELAX_screen_hole_depth(bare, {'E1'}, 1.4);
    assert(isnan(hole.depth) && hole.note == "Fewer than two electrodes have coordinates");

    hole = RELAX_screen_hole_depth(struct('labels', {'E1', 'E2'}), {'E1'}, 1.4);
    assert(isnan(hole.depth) && hole.note == "Fewer than two electrodes have coordinates");
    fprintf('test_RELAX_screen_hole_depth: PASS\n');
end

function locs = lattice()
    [columns, rows] = meshgrid(-3:3, -3:3);
    isCentre = columns(:) == 0 & rows(:) == 0;
    [~, centreFirst] = sort(~isCentre);
    columns = columns(centreFirst);
    rows = rows(centreFirst);
    labels = compose('E%d', (1:numel(columns))');
    x = columns + 0.5 * mod(rows, 2);
    y = rows * sqrt(3) / 2;
    locs = struct('labels', labels, 'X', num2cell(x), 'Y', num2cell(y), 'Z', 0);
end

function labels = labelsWithin(locs, radius)
    centre = locs(strcmp({locs.labels}, 'E1'));
    distance = hypot([locs.X] - centre.X, [locs.Y] - centre.Y);
    labels = {locs(distance <= radius + 1e-9).labels};
end

function labels = rowOfFive(locs)
    onRow = abs([locs.Y]) < 1e-9 & abs([locs.X]) <= 2;
    labels = {locs(onRow).labels};
end
