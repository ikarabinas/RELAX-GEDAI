% RELAX_screen_write_csv() - Write the screen table as a CSV, with logical columns as 0/1.
%
% Usage:
%   >> RELAX_screen_write_csv(screen, csvPath);

function RELAX_screen_write_csv(screen, csvPath)
    forCsv = screen;
    for column = string(forCsv.Properties.VariableNames)
        if islogical(forCsv.(column)), forCsv.(column) = double(forCsv.(column)); end
    end
    writetable(forCsv, csvPath);
end
