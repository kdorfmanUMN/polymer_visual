% This function accepts arrays R, x, y, and z representing the data
% and gridpoint coordinates for an SCFT solution in a single unit
% cell, and generates new arrays R2, x2, y2, and z2 with different
% boundaries as specified by inputs alim, blim, and clim. Periodic
% boundary conditions are assumed in order to extend outside of a 
% single unit cell. If the cell limit defined by the user is not on 
% a gridpoint, we round that cell limit to the nearest value that 
% would land on a gridpoint.

function [R2,x2,y2,z2] = change_cell_lims(R,x,y,z,options)

    arguments
        R
        x
        y
        z
        options.alim = [0,1];
        options.blim = [0,1];
        options.clim = [0,1];

        % normalVec is used only when plotting thin films, so that we 
        % don't enforce periodic behavior in the direction normal to
        % the thin film walls/substrate. If x is normal to the walls
        % then normalVec = 0, if y is normal to the walls then 
        % normalVec = 1, and if z is normal to the walls then normalVec
        % = 2. A value of -1 (the default) indicates that there is not
        % a thin film boundary, and we assume 3D periodicity.
        % normalVec indexes from 0 rather than 1 to match the input file
        % structure of pscfpp.
        options.normalVec = -1;

    end

    % If alim, blim, and clim are all [0,1], no modifications are needed
    if isequal(options.alim,[0,1]) && isequal(options.blim,[0,1]) && ...
       isequal(options.clim,[0,1])
        R2 = R; x2 = x; y2 = y; z2 = z;
        return
    end

    normalVec = options.normalVec;
    
    % Get basis vectors and grid array
    basis = [x(end,1,1)-x(1,1,1), y(end,1,1)-y(1,1,1), z(end,1,1)-z(1,1,1);
             x(1,end,1)-x(1,1,1), y(1,end,1)-y(1,1,1), z(1,end,1)-z(1,1,1);
             x(1,1,end)-x(1,1,1), y(1,1,end)-y(1,1,1), z(1,1,end)-z(1,1,1)];
    grid_bulk = size(x)-1;
    grid = grid_bulk;

    % determine alim, blim, and clim in grid coordinates
    alim_grid = round(options.alim * grid(1));
    blim_grid = round(options.blim * grid(2));
    clim_grid = round(options.clim * grid(3));

    % Check if this is a thin film and correct grid if necessary
    if normalVec == 0
        grid(1) = grid(1) + 1;
    elseif normalVec == 1
        grid(2) = grid(2) + 1;
    elseif normalVec == 2
        grid(3) = grid(3) + 1;
    end
    
    % Create empty arrays to store new data
    R2 = zeros(alim_grid(2)-alim_grid(1)+1, blim_grid(2)-blim_grid(1)+1,...
               clim_grid(2)-clim_grid(1)+1, size(R,4));
    x2 = zeros(alim_grid(2)-alim_grid(1)+1, blim_grid(2)-blim_grid(1)+1,...
               clim_grid(2)-clim_grid(1)+1);
    y2 = zeros(alim_grid(2)-alim_grid(1)+1, blim_grid(2)-blim_grid(1)+1,...
               clim_grid(2)-clim_grid(1)+1);
    z2 = zeros(alim_grid(2)-alim_grid(1)+1, blim_grid(2)-blim_grid(1)+1,...
               clim_grid(2)-clim_grid(1)+1);
    
    % Loop through all gridpoints in the new cell limits and
    % store correct data in R2, x2, y2, and z2.
    ix = 0; iy = 0; iz = 0;
    for xg = alim_grid(1):alim_grid(2)
        ix = ix + 1;
        x_cell = fit_in_cell(xg+1,grid(1));
        for yg = blim_grid(1):blim_grid(2)
            iy = iy + 1;
            y_cell = fit_in_cell(yg+1,grid(2));
            for zg = clim_grid(1):clim_grid(2)
                iz = iz + 1;
                z_cell = fit_in_cell(zg+1,grid(3));
    
                R2(ix,iy,iz,:) = R(x_cell,y_cell,z_cell,:);
                
                xval = (basis(1,1) * xg/grid_bulk(1)) + ...
                       (basis(2,1) * yg/grid_bulk(2)) + ...
                       (basis(3,1) * zg/grid_bulk(3)) + x(1,1,1);
                yval = (basis(1,2) * xg/grid_bulk(1)) + ...
                       (basis(2,2) * yg/grid_bulk(2)) + ...
                       (basis(3,2) * zg/grid_bulk(3)) + x(1,1,1);
                zval = (basis(1,3) * xg/grid_bulk(1)) + ...
                       (basis(2,3) * yg/grid_bulk(2)) + ...
                       (basis(3,3) * zg/grid_bulk(3)) + x(1,1,1);

                % If thin film, leave gap where wall should be
                if normalVec == 0
                    xval = xval + (floor(xg/grid(1)) * 2 * x(1,1,1));
                    xval = xval - (floor(xg/grid(1))*basis(1,1)/grid(1));
                elseif normalVec == 1
                    yval = yval + (floor(yg/grid(2)) * 2 * y(1,1,1));
                    yval = yval - (floor(yg/grid(2))*basis(2,2)/grid(2));
                elseif normalVec == 2
                    zval = zval + (floor(zg/grid(3)) * 2 * z(1,1,1));
                    zval = zval - (floor(zg/grid(3))*basis(3,3)/grid(3));
                end
    
                x2(ix,iy,iz) = round(xval,10);
                y2(ix,iy,iz) = round(yval,10);
                z2(ix,iy,iz) = round(zval,10);
            end
            iz = 0;
        end
        iy = 0;
    end

end