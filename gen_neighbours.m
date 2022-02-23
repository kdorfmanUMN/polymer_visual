function neighbours = gen_neighbours(cell_d,angles,phase)
    basis = get_basis(cell_d,angles);
    atomidx = get_atomloc(phase);
    neighbour_idx = [1,0,0;0,1,0;0,0,1;1,1,0;1,-1,0;1,0,1;1,0,-1;0,1,1;0,1,-1;1,1,1;-1,1,1;1,-1,1;1,1,-1];
    
    neighboursTemp = atomidx;
    for i = 1:size(neighbour_idx,1)
        neighboursTemp = cat(1,neighboursTemp, atomidx + neighbour_idx(i,:));
        neighboursTemp = cat(1,neighboursTemp, atomidx - neighbour_idx(i,:));    
    end
    neighbours = neighboursTemp; % if you just want to return the idx
    %neighbours = round(((basis')*(neighboursTemp'))',7);
end