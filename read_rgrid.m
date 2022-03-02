function [R,x,y,z,dim,lattype,cell_d,angles,n_mnr,grid] = read_rgrid(filename)

    tmp = fopen(filename);
    C = textscan(tmp,'%s','delimiter', '\n');
    C=C{1};
    fclose(tmp); clear tmp;

    ic = 0;
    start_row = 1000;
    % Reading the Preamble of the File
    
    while ic <= start_row
        ic = ic +1;
        
        % This section is a way of checking when the data actually starts...
        % checks if a line, and the line below it, sum to 1. Does this 5 times.
        % If so, data has probably started.
        
        % interesting.
        
        
        %if round(sum (sscanf(C{ic},'%f')),2)== 1.00 && round(sum (sscanf(C{ic+1},'%f')),2)== 1.00
        %    ndata = ndata+1;
        %else
        %    ndata = 0;
        %end
        
        % this section reads in preamble.
        
        if strcmp(strrep(char(C(ic)),' ', ''),'dim')==1
            dim = str2double(C{ic+1});          % Reads the grid dimensions
        elseif strcmp(strrep(char(C(ic)),' ', ''),'crystal_system')==1
            lattype = strrep(C(ic+1), '''', '');   % Reads the system type
        elseif strcmp(strrep(char(C(ic)),' ', ''),'cell_param')==1
            param = sscanf(C{ic+1},'%f')';            % Reads the cell parameters
        elseif strcmp(strrep(char(C(ic)),' ', ''),'N_monomer')==1
            n_mnr = str2double(C{ic+1});        % Reads the number of monomers
        elseif strcmp(strrep(char(C(ic)),' ', ''),'ngrid')==1
            grid = sscanf(C{ic+1},'%f')';            % Reads the grid size
            start_row = ic+2;       % Records the row in which the volume fractions start
        end
    end
    
    end_info = start_row - 1;                % Records the row in which the supplementary information ends
                     
    
    % Reading the Grid Points From the File, but not in grid fashion yet..
    
    A = zeros(length(C) - end_info,n_mnr);
    
    for i =start_row:length(C)
        A(i - end_info,:) = sscanf(C{i},'%f')';
    end
    
    
    % Get x,y,z grid points and unit cell parameters
    
    [x,y,z,cell_d,angles] = gen_xyz(lattype,param,grid);
    
    % place points from A on grid, for each monomer. Finally. R is a 4D array.
    % 3D are for the x,y,z grid. The 4th D is for the monomer type.
    
    R = rearrangePts(A, grid, dim);

end

function R = rearrangePts(A, grid, dim)

    n_mnr = size(A,2);
    R = zeros([grid+1 n_mnr]); % store values of volume fraction on grid
    counter = 0;
    
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

function [x,y,z,cell_d,angle] = gen_xyz(lattype, param, grid)

    % get crystal system dimensions and angles
    
    if strcmp(lattype,'hexagonal') == 1
        angle = [pi/2 pi/2 (2*pi)/3];
        cell_d = param;
    elseif strcmp(lattype,'cubic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'tetragonal') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'orthorhombic') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    elseif strcmp(lattype,'triclinic') == 1
        angle = [param(4) param(5) param(6)];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'monoclinic') == 1
        angle = [pi/2 param(4) pi/2];
        cell_d = [param(1) param(2) param(3)];
    elseif strcmp(lattype,'trigonal') == 1
        angle = [param(2) param(2) param(2)];
        cell_d = param(1);
    elseif strcmp(lattype,'lamellar') == 1
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    else
        angle = [pi/2 pi/2 pi/2];
        cell_d = param;
    end
    
    % extract unit cell dimensions, such that there are 3 for each system
    if(length(cell_d)==1)
        new_cell(1:3) = cell_d;             % Cubic crystals
    elseif(length(cell_d)==2)
        new_cell(1:2) = cell_d(1);            % Tetragonal crystals
        new_cell(3)   = cell_d(2);
    else
        new_cell = cell_d;                    % Orthorhombic crystals
    end

    clear cell_d; cell_d = new_cell;

    if(length(grid)==1)
        grid(2) = grid(1);                   % 3D grid for 1D crystals
        grid(3) = grid(1);
    elseif(length(grid)==2)
        grid(3) = grid(1);                  % 3D grid for 2D crystals
    end
    
    basis = get_basis(cell_d,angle);

    % matrices for grid coords
    x = zeros(grid); % x coords on grid
    y = zeros(grid); % y coords on grid
    z = zeros(grid); % z coords on grid
    
    nround = 10;
    
    % calculate the coordinates in x,y,z.
    for iz=1:grid(3)+1
        for iy=1:grid(2)+1
            for ix=1:grid(1)+1
                xtemp = (basis(1,1) * (ix-1)/grid(1)) + (basis(2,1) * (iy-1)/grid(2)) + (basis(3,1) * (iz-1)/grid(3));
                ytemp = (basis(2,2) * (iy-1)/grid(2)) + (basis(3,2) * (iz-1)/grid(3));
                ztemp = (basis(3,3) * (iz-1)/grid(3));

                x(ix,iy,iz) = round(xtemp,nround);
                y(ix,iy,iz) = round(ytemp,nround);
                z(ix,iy,iz) = round(ztemp,nround);
            end
        end
    end
end