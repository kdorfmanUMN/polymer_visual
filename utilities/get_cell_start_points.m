% The purpose of this function is to take an lower and upper bound of a
% discretized 1D unit cell in reduced coordinates (e.g., [-0.5,2]), as well
% as an integer indicating the number of gridpoints in the domain [0,1),
% and return an array containing the index of the gridpoint at the start of
% each new unit cell in the domain, as well as the first and last 
% gridpoint. This is easier to understand with an example and an
% explanation of how/why we use it.

% The simplest use case for this is if we have an array containing multiple
% unit cells of periodic data, and we want to plot each unit cell
% separately (perhaps to make each one a distinct color). We know the array
% spans the domain [-0.5,2] in reduced coordinates, and each unit cell has
% 64 gridpoints. In order to plot each one separately, we need to know the
% index of the starting point of each partial or full unit cell in the 
% domain, which in this case would be [1,33,97,161] corresponding to the 
% reduced coordinate points [-0.5,0,1,2]. The array [1,33,97,161] is what
% would be output by this function. 

% The lower and upper bound in reduced coordinates is input as a
% two-element array "lims", and the number of gridpoints per unit cell is
% input as an integer "nsteps".

function out = get_cell_start_points(lims,nsteps)
    
    arguments
        lims   (1,2) {mustBeNumeric} % Lower/upper bound in reduced coords
        nsteps (1,1) {mustBeNumeric} % # of gridpoints per unit cell
    end
    
    lims = round(lims,8); 
    out = floor(lims(1)):ceil(lims(2));
    out(1) = lims(1); out(end) = lims(2);
    out = out - lims(1);
    out = out * nsteps + 1;

end