% This function accepts arrays R, x, y, and z representing the data
% and gridpoint coordinates for an SCFT solution in a single unit
% cell, and generates new arrays R2, x2, y2, and z2 with different
% boundaries as specified by inputs xlim, ylim, and zlim. Periodic
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
        options.xlim = [0,1];
        options.ylim = [0,1];
        options.zlim = [0,1];

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

    % If xlim, ylim, and zlim are all [0,1], no modifications are needed
    if isequal(options.xlim,[0,1]) && isequal(options.ylim,[0,1]) && ...
       isequal(options.zlim,[0,1])
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

    % Check if this is a thin film and correct grid if necessary
    if normalVec == 0
        grid(1) = grid(1) + 1;
    elseif normalVec == 1
        grid(2) = grid(2) + 1;
    elseif normalVec == 2
        grid(3) = grid(3) + 1;
    end
    
    % determine xlim, ylim, and zlim in grid coordinates
    xlim_grid = round(options.xlim * grid(1));
    ylim_grid = round(options.ylim * grid(2));
    zlim_grid = round(options.zlim * grid(3));
    
    % Create empty arrays to store new data
    R2 = zeros(xlim_grid(2)-xlim_grid(1)+1, ylim_grid(2)-ylim_grid(1)+1,...
               zlim_grid(2)-zlim_grid(1)+1, size(R,4));
    x2 = zeros(xlim_grid(2)-xlim_grid(1)+1, ylim_grid(2)-ylim_grid(1)+1,...
               zlim_grid(2)-zlim_grid(1)+1);
    y2 = zeros(xlim_grid(2)-xlim_grid(1)+1, ylim_grid(2)-ylim_grid(1)+1,...
               zlim_grid(2)-zlim_grid(1)+1);
    z2 = zeros(xlim_grid(2)-xlim_grid(1)+1, ylim_grid(2)-ylim_grid(1)+1,...
               zlim_grid(2)-zlim_grid(1)+1);
    
    % Loop through all gridpoints in the new cell limits and
    % store correct data in R2, x2, y2, and z2.
    ix = 0; iy = 0; iz = 0;
    for xg = xlim_grid(1):xlim_grid(2)
        ix = ix + 1;
        x_cell = fit_in_cell(xg+1,grid(1));
        for yg = ylim_grid(1):ylim_grid(2)
            iy = iy + 1;
            y_cell = fit_in_cell(yg+1,grid(2));
            for zg = zlim_grid(1):zlim_grid(2)
                iz = iz + 1;
                z_cell = fit_in_cell(zg+1,grid(3));
    
                R2(ix,iy,iz,:) = R(x_cell,y_cell,z_cell,:);
                
                xval = (basis(1,1) * xg/grid_bulk(1)) + ...
                       (basis(2,1) * yg/grid_bulk(2)) + ...
                       (basis(3,1) * zg/grid_bulk(3)) + x(1,1,1);
                yval = (basis(2,2) * yg/grid_bulk(2)) + ...
                       (basis(3,2) * zg/grid_bulk(3)) + y(1,1,1);
                zval = (basis(3,3) * zg/grid_bulk(3)) + z(1,1,1);

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