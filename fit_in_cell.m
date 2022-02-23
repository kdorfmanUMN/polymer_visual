% exists to correct for any grid coordinate that is outside the bounds of a
% unit cell. Accepts any grid coordinate, and returns the grid coordinate
% *inside* the unit cell that matches that grid coordinate. For example, if
% the grid is discretized into 64 steps in a given Cartesian direction
% (gridmax = 64), and we have a coordinate at -30 (val = -30), this
% function will output a value of 34. If the coordinate is 288, the
% function will output 32.

% Note: this function is for plotting, which means that we include data on
% all faces of the unit cell. If fixed_val = gridmax + 1, this means that
% fixed_val is on the far face of the unit cell, which we consider as being
% "in the unit cell" in this function.

% Also note: We are using matlab indexing here, starting from 1. So, if
% fixed_val = 1, that's in the unit cell; if fixed_val = 0, it's outside.

function fixed_val = fit_in_cell(val,gridmax)
    fixed_val = val;
    if fixed_val < 1
        while fixed_val < 1
            fixed_val = fixed_val + gridmax;
        end
    elseif fixed_val > gridmax + 1
        while fixed_val > gridmax + 1
            fixed_val = fixed_val - gridmax;
        end
    end
end