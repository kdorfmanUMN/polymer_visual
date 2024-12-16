%% Main function: read_rgrid

% Reads data from a PSCF r-grid file (filename) and returns arrays R, x, y,
% and z that contain the necessary data to plot the composition profiles. R
% contains the composition of each species at each gridpoint, so R(i,j,k,l)
% is the composition of species l at point (i,j,k). x, y, and z contain the
% x, y and z coordinates of each gridpoint, respectively, so x(i,j,k)
% contains the x-coordinate of gridpoint (i,j,k). Data in 1D or 2D are
% extended to 3D. Each output array contains data for one additional
% gridpoint in each dimension, to show the "back side" of the periodic unit
% cell as identical to the "front side" that contacts gridpoint (1,1,1).

% The function also returns dim, an integer (either 1, 2, or 3)
% representing the dimensionality of the data in the original r-grid file.
% It further returns lattype, a string representing the crystal system of
% the lattice (e.g., "orthorhombic").

% The r-grid file of an FTS trajectory will, in general, contain multiple
% fields, labeled by an index i in the field file that begins at 0 and
% increases with each consecutive field. If an FTS trajectory is provided
% as the input file to this function, the field that is read will, by
% default, be the first field in the file, with i = 0. An optional input
% parameter fieldId may be provided, in which case the field with index
% fieldId will be read instead. 

function [R, x, y, z, dim, lattype] = read_rgrid(filename, fieldId)

    arguments

        % String that represents the path to the rgrid file containing the 
        % data that will be read
        filename;

        % Optional index that indicates which field to read from an FTS
        % trajectory. If the input file is not an FTS trajectory, this
        % value is not used. Default is 0 (reads the first field in the
        % file).
        fieldId = 0;

    end
    
    % Ensure that the code below can access our utilities
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(filepath+"/utilities")

    % Open file, store data
    tmp = fopen(filename);
    C = textscan(tmp,'%s','delimiter', '\n');
    C=C{1};
    fclose(tmp); clear tmp;
        
    % Check that field file header is in the expected format
    if (~strcmp(strip(C{2}), 'dim') || ...
        ~strcmp(strip(C{4}), 'crystal_system') || ...
        ~strcmp(strip(C{8}), 'cell_param'))
        error("Field file header cannot be parsed.");
    end

    % Read necessary data from field file header
    dim = str2double(C{3}); 
    lattype = strrep(C(5), '''', ''); 
    param = sscanf(C{9},'%f')';

    % Find N_monomer
    sym = true;
    if strcmp(strip(C{10}), 'N_monomer') % no symmetry
        n_mnr = str2double(C{11}); 
        sym = false; 
    elseif strcmp(strip(C{12}), 'N_monomer') % with symmetry
        n_mnr = str2double(C{13}); 
    else
        error("N_monomer not found.")
    end

    % Determine if this is an FTS simulation field file
    fts = false;
    if strcmp(C{13}, 'i = 0')
        fts = true;
    end

    j = 13;
    if fts
        % If this is an FTS simulation, find field with index fieldId
        fieldIdHeader = sprintf('i = %d',fieldId);
        lenC = length(C);
        while ~strcmp(C{j},fieldIdHeader)
            j = j + 1;
            if j > lenC
                error("Field with index %d not found.",fieldId);
            end
        end
        j = j + 1;
    elseif sym
        j = 14;
    else 
        j = 12;
    end

    if strcmp(strrep(char(C(j)),' ', ''),'ngrid') || ...
       strcmp(strrep(char(C(j)),' ', ''),'mesh')
        % Reads the grid size
        grid = sscanf(C{j+1},'%f')'; 
        start_row = j+2; % the row where volume fractions start
    else
        error("Mesh dimensions not found");
    end
        
    % Extend 1D and 2D grids to a 3D grid
    if(isscalar(grid)) % grid is only one element
        grid(2) = 3; %grid(1);
        grid(3) = 3; %grid(1);
    elseif(length(grid) == 2)
        grid(3) = 3; %grid(1);
    end
    
    % Read grid points from the file into a linear array
    n_pts = prod(grid);
    A = zeros(n_pts,n_mnr);
    
    for i =1:n_pts
        A(i,:) = sscanf(C{i+start_row-1},'%f')';
    end
    
    % Get x,y,z grid points and unit cell parameters
    [x,y,z,~] = gen_xyz(lattype,param,grid);
    
    % place points from A on grid, for each monomer. R is a 4D array.
    % 3D are for the x,y,z grid. The 4th D is for the monomer type.
    R = rearrange_pts(A, grid, dim);

end

%% rearrange_pts

% This function takes a linear array A with dimensions (n,n_mnr) (where n
% is the number of gridpoints and n_mnr is the number of data values at
% each gridpoint) and rearranges it into a 4D array R where the data are
% indexable by gridpoint, e.g. R(x,y,z,:) is the data at gridpoint (x,y,z).
% The underlying grid is specified by the "grid" input (e.g. [32 32 32]),
% and the dimension of the incoming data is specified by "dim." If the data
% are in 1D or 2D and higher dimensions are specified in "grid" (e.g. the
% data are in 1D with 32 gridpoints but grid=[32 32 32]) then the data will
% be duplicated along all higher dimensions to create an output array of
% the same dimensionality as "grid".

% Additionally, this function extends the data by one gridpoint in each
% dimension, so a 32x32x32 grid will generate an R array of dimensions
% 33x33x33xn_mnr. The data in these additional gridpoints are determined by
% assuming periodic boundary conditions. This extension of the data allows
% us to visualize the unit cell and see the periodic behavior along each
% dimension clearly, since this extension results in all opposite faces of
% the unit cell being identical.
function R = rearrange_pts(A, grid, dim)

    n_mnr = size(A,2);
    R = zeros([grid+1 n_mnr]); % store values of volume fraction on grid
    counter = 0;

    if length(grid) == 1
        grid(2) = 1; grid(3) = 1;
    elseif length(grid) == 2
        grid(3) = 1;
    elseif length(grid) > 3
        error("grid vector can not contain more than 3 values")
    end
    
    for iz=1:grid(3)+1
        for iy=1:grid(2)+1
            for ix=1:grid(1)+1
                counter = counter + 1;
                for in = 1:n_mnr
                    if ix == grid(1)+1
                        R(grid(1)+1,:,:,in) = R(1,:,:,in); % periodic bc's
                        counter = counter - (1/n_mnr); 
                    elseif iy == grid(2) + 1
                        R(:,grid(2)+1,:,in) = R(:,1,:,in);
                        counter = counter - (1/n_mnr);
                    elseif iz == grid(3) + 1
                        R(:,:,grid(3)+1,in) = R(:,:,1,in);
                        counter = counter - (1/n_mnr);
                    else
                        R(ix,iy,iz,in) = A(round(counter),in);
                    end
                end
            end
            if(dim==1)
                counter = 0;
            end
        end
        if(dim==2)
            counter=0;
        end
    end
end

%% gen_xyz:

% This function uses the crystal system (lattype) and the lattice 
% parameters (param) to determine all lattice vector lengths (cell_d) 
% and unit cell angles (angle). This is used to generate a 3x3 array
% containing the lattice basis vectors (see get_basis function below).
% Finally, these lattice basis vectors are used to generate arrays x, 
% y, and z containing the coordinates of each gridpoint, where the 3D 
% mesh (e.g. 32x32x32) is input as a 3-element vector (grid). Generalized
% for any crystal system.

function [x,y,z,basis] = gen_xyz(lattype, param, grid)

    % get crystal system dimensions and angles
    
    if strcmp(lattype,'hexagonal') == 1
        angle = [pi/2 pi/2 (2*pi)/3];
        if length(param) == 2 % assume 3d
            cell_d = [param(1) param(1) param(2)];
        else % assume 2d, param(2) does not exist
            cell_d = [param(1) param(1) param(1)];
        end
    elseif strcmp(lattype,'cubic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param param param];
    elseif strcmp(lattype,'tetragonal') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param(1) param(1) param(2)];
    elseif strcmp(lattype,'orthorhombic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'triclinic') == 1
        phi = param(4);
        theta = param(5);
        gamma = param(6);
        
        % convert phi and theta to alpha and beta:
        beta = acos(sin(theta)*cos(phi));
        alpha = acos(sin(theta)*sin(phi) - ...
                     (cos(beta)*cos(gamma)/sin(gamma)));

        angle = [alpha beta gamma];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'monoclinic') == 1
        angle = [pi/2 param(4) pi/2];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'trigonal') == 1
        angle = [param(2) param(2) param(2)];
        cell_d = [param(1) param(1) param(1)];
    elseif strcmp(lattype,'rhombohedral') == 1
        angle = [param(2) param(2) param(2)];
        cell_d = [param(1) param(1) param(1)];
    elseif strcmp(lattype,'square') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param(1) param(1) param(1)];
    elseif strcmp(lattype,'rectangular') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param(1) param(2) param(1)];
    elseif strcmp(lattype,'oblique') == 1
        angle = [pi/2 pi/2 param(3)];
        cell_d = [param(1) param(2) param(1)];
    elseif strcmp(lattype,'rhombic') == 1
        angle = [pi/2 pi/2 param(2)];
        cell_d = [param(1) param(1) param(1)];
    elseif strcmp(lattype,'lamellar') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = [param(1) param(1) param(1)];
    else
        error("unknown crystal system")
    end

    % Extend 1D and 2D grids to a 3D grid
    if(length(grid)==1)
        grid(2) = grid(1);
        grid(3) = grid(1);
    elseif(length(grid)==2)
        grid(3) = grid(1);
    end
    
    % Get lattice basis vectors
    basis = get_basis(cell_d,angle);

    % matrices for grid coords
    x = zeros(grid+1); % x coords on grid
    y = zeros(grid+1); % y coords on grid
    z = zeros(grid+1); % z coords on grid
    
    nround = 10;
    
    % calculate the coordinates in x,y,z.
    for iz=1:grid(3)+1
        for iy=1:grid(2)+1
            for ix=1:grid(1)+1
                xtemp = (basis(1,1) * (ix-1)/grid(1)) + ...
                        (basis(2,1) * (iy-1)/grid(2)) + ...
                        (basis(3,1) * (iz-1)/grid(3));
                ytemp = (basis(2,2) * (iy-1)/grid(2)) + ...
                        (basis(3,2) * (iz-1)/grid(3));
                ztemp = (basis(3,3) * (iz-1)/grid(3));

                x(ix,iy,iz) = round(xtemp,nround);
                y(ix,iy,iz) = round(ytemp,nround);
                z(ix,iy,iz) = round(ztemp,nround);
            end
        end
    end
end
