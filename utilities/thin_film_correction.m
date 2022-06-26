% This function applies a correction to the composition data of an SCFT
% solution that was calculated under a thin film constraint, so that the
% sum of the compositions do not decay to zero when moving toward a wall,
% but rather stays fixed at 1. The function also deletes any data points
% that are "in the wall" (anywhere where the volume fraction of the wall is
% >0.5). The edited arrays R,x,y,z can then be plotted and it will look
% like a polymer system confined between walls, rather than having a small
% decay to 0 concentration for each species at each side.

% The inputs normalVec, T, and t 

function [R,x,y,z,basis] = thin_film_correction(R,x,y,z,basis,normalVec,...
                                                t,T,rotate)

    % Setup
    grid = size(x);
    if ~exist('rotate','var')
        rotate = true; % default behavior for rotate
    end
    delete = [];

    % Rotate data if requested, so that the wall is normal to the z axis
    if rotate
        xf = reshape(x,[],1);
        yf = reshape(y,[],1);
        zf = reshape(z,[],1);
        
        % get rotation matrix
        if normalVec == 0
            rot = [0 0 1; 0 1 0; -1 0 0]; % rotate 90° around y axis
        elseif normalVec == 1
            rot = [1 0 0; 0 0 -1; 0 1 0]; % rotate 90° around x axis
        else % normalVec == 2
            rot = [1 0 0; 0 1 0; 0 0 1];  % do nothing, already rotated
        end

        new_coords = (rot * ([xf,yf,zf]'))'; % apply rotation matrix
        basis = (rot * basis')'; % rotate basis vectors too
        
        % store rotated coords in x, y, and z matrices
        x = reshape(new_coords(:,1),size(x));
        y = reshape(new_coords(:,2),size(y));
        z = reshape(new_coords(:,3),size(z));
    end

    L = norm(basis(normalVec+1,:));

    % Correct data so that the sum of all polymer species adds to 1
    for ix = 1:grid(1)
        for iy = 1:grid(2)
            for iz = 1:grid(3)
                pos = [ix,iy,iz];
                d = (pos(normalVec+1)-1) * L / (grid(normalVec+1)-1);
                rho_w = 0.5*(1+tanh(4*(((.5*(T-L))+abs(d-(L/2)))/t)));
                
                if rho_w <= 0.5
                    R(ix,iy,iz,:) = R(ix,iy,iz,:) ./ (1-rho_w);
                elseif ~ismember(pos(normalVec+1),delete)
                    delete(end+1) = pos(normalVec+1); %#ok<AGROW> 
                end
                
            end
        end
    end
    
    % Delete all data that is "in the wall"
    if normalVec == 0
        R(delete,:,:,:) = [];
        x(delete,:,:) = [];
        y(delete,:,:) = [];
        z(delete,:,:) = [];
    elseif normalVec == 1
        R(:,delete,:,:) = [];
        x(:,delete,:) = [];
        y(:,delete,:) = [];
        z(:,delete,:) = [];
    else % normalVec == 2
        R(:,:,delete,:) = [];
        x(:,:,delete) = [];
        y(:,:,delete) = [];
        z(:,:,delete) = [];
    end
    
end