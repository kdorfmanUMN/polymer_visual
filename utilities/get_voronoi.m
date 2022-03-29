function [v,c] = get_voronoi(basis,phase)
    % get number of atoms 
    atomidx = get_atomloc(phase);
    n_atoms = size(atomidx,1);
    % generate large number of neighbours
    neighbours = gen_neighbours(basis,phase);
    % find voronoi cells for all points
    [v,c] = voronoin(neighbours); % v is 
    % reduce to only the cells relevant to the unit cell points, not their
    % neighbours
    %c = c(1:n_atoms);
end