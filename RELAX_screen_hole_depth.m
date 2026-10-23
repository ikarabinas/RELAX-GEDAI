% RELAX_screen_hole_depth() - How deep the deepest hole of interpolated electrodes is.
%
% Usage:
%   >> hole = RELAX_screen_hole_depth(chanlocs, interpolatedLabels, limit);
%
% Inputs:
%   chanlocs           - EEG.allchan-style struct array with labels, X, Y and Z; any coordinate may be []
%   interpolatedLabels - cell array or string array of the interpolated electrode labels
%   limit              - an interpolated electrode is listed in hole.electrodes when its nearest real
%                        electrode is further than limit * spacing
%
% Outputs:
%   hole.depth      - distance from the interpolated electrode furthest from any real electrode to its
%                     nearest real electrode, in spacings. 0 when nothing is interpolated, NaN when it
%                     cannot be measured
%   hole.deepest    - label of that electrode, "" when nothing is interpolated or depth is NaN
%   hole.electrodes - [1 x k] string, the interpolated electrodes beyond limit, in chanlocs order
%   hole.spacing    - median over electrodes with coordinates of the distance to the nearest other
%                     electrode, in coordinate units (NaN with fewer than two such electrodes)
%   hole.note       - why depth is NaN, "" otherwise

function hole = RELAX_screen_hole_depth(chanlocs, interpolatedLabels, limit)

    labels = string({chanlocs.labels});
    isInterpolated = ismember(labels, string(interpolatedLabels));
    coordinates = [coordinateOf(chanlocs, 'X'), coordinateOf(chanlocs, 'Y'), coordinateOf(chanlocs, 'Z')];
    hasCoordinates = all(isfinite(coordinates), 2)';

    hole = struct('depth', 0, 'deepest', "", 'electrodes', strings(1, 0), 'spacing', NaN, 'note', "");
    hole.spacing = electrodeSpacing(coordinates(hasCoordinates, :));

    if ~any(isInterpolated)
        return
    end
    hole.depth = NaN;

    realWithCoordinates = ~isInterpolated & hasCoordinates;
    interpolatedWithCoordinates = isInterpolated & hasCoordinates;
    if nnz(hasCoordinates) < 2
        hole.note = "Fewer than two electrodes have coordinates";
        return
    elseif ~any(realWithCoordinates)
        hole.note = "No real electrode has coordinates";
        return
    elseif any(isInterpolated & ~hasCoordinates)
        hole.note = "Interpolated electrodes without coordinates: " + strjoin(labels(isInterpolated & ~hasCoordinates), ", ");
        return
    end

    nearestRealDistance = NaN(size(labels));
    nearestRealDistance(interpolatedWithCoordinates) = min(pdist2(coordinates(interpolatedWithCoordinates, :), coordinates(realWithCoordinates, :)), [], 2);
    [deepestDistance, deepestIndex] = max(nearestRealDistance);
    hole.depth = deepestDistance / hole.spacing;
    hole.deepest = labels(deepestIndex);
    hole.electrodes = labels(nearestRealDistance > limit * hole.spacing);
end

function values = coordinateOf(chanlocs, field)
    values = NaN(numel(chanlocs), 1);
    if ~isfield(chanlocs, field)
        return
    end
    present = ~cellfun(@isempty, {chanlocs.(field)});
    values(present) = [chanlocs(present).(field)];
end

function spacing = electrodeSpacing(coordinates)
    if size(coordinates, 1) < 2
        spacing = NaN;
        return
    end
    distances = pdist2(coordinates, coordinates);
    distances(logical(eye(size(distances)))) = Inf;
    spacing = median(min(distances, [], 2));
end
