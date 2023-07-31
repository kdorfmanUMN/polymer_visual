% This function applies a correction to the composition data of an SCFT
% solution that was calculated under a thin film constraint, so that the
% sum of the compositions do not decay to zero when moving toward a wall,
% but rather stays fixed at 1. The function also deletes any data points
% that are "in the wall" (anywhere where the volume fraction of the wall is
% >0.5). The edited arrays R,x,y,z can then be plotted and it will look
% like a polymer system confined between walls, rather than having a small
% decay to 0 concentration for each species at each side.

% Aside from the original R,x,y,z arrays, the basis vector array (output by
% utilities/get_basis.m) is also a required input, as are the thin film
% parameters normalVec, t, T corresponding to normalVec,
% interfaceThickness, and wallThickness in the param file for a thin film.

% An optional parameter, rotate, is a boolean which indicates whether
% the data should be adjusted so that the wall is in the x-y plane. If
% rotate=true, this adjustment will be applied. Otherwise, R,x,y,z arrays
% will not be rotated. (note: technically the transformation we apply is a 
% permutation, not a rotation, but we call it a rotation anyway).

% Finally, an optional input keep_nan is a boolean. If true, the output R
% will be the same shape as the original R array, but with NaN as the value
% at every point that is "in the wall". Otherwise, if this is false (the
% default), the data "in the wall" will instead just be deleted from R, x,
% y, and z. keep_nan=true is used in line_profile and contour_plot.

function [R,x,y,z,basis] = thin_film_correction(R,x,y,z,normalVec,t,T,...
                                                rotate,keep_nan)
    arguments
        % Required parameters
        R               {mustBeNumeric} % Data
        x               {mustBeNumeric} % x-coordinates
        y               {mustBeNumeric} % y-coordinates
        z               {mustBeNumeric} % z-coordinates
        normalVec (1,1) {mustBeNumeric} % Which lattice vector is normal to
                                        % the wall? 0 (a), 1 (b), or 2 (c)
        t         (1,1) {mustBeNumeric} % interface thickness
        T         (1,1) {mustBeNumeric} % wall thickness
        
        % Optional parameters
        rotate   (1,1) {mustBeNumericOrLogical} = true  % rotate unit cell?
        keep_nan (1,1) {mustBeNumericOrLogical} = false % keep NaN values
                                                        % (T) or trim R, x, 
                                                        % y, z arrays (F)?
    end

    % Setup
    basis = [x(end,1,1),y(end,1,1),z(end,1,1);
             x(1,end,1),y(1,end,1),z(1,end,1);
             x(1,1,end),y(1,1,end),z(1,1,end)];
    if ~exist('rotate','var')
        rotate = true; % default behavior for rotate
    end

    delete = [];

    % Rotate data if requested, so that the wall is normal to the z axis
    % (actually, we are permuting the coords, but it's a similar effect)
    if rotate
        if normalVec == 0
            R = permute(R,[2,3,1,4]);
            tmp = x;
            x = permute(y,[2,3,1]);
            y = permute(z,[2,3,1]);
            z = permute(tmp,[2,3,1]);
            clear tmp;
            basis = basis([2,3,1],[2,3,1]);
        elseif normalVec == 1
            R = permute(R,[3,1,2,4]);
            tmp = x;
            x = permute(z,[3,1,2]);
            z = permute(y,[3,1,2]);
            y = permute(tmp,[3,1,2]);
            clear tmp;
            basis = basis([3,1,2],[3,1,2]);
        end % if normalVec == 2 do nothing
        normalVec = 2;
    end
    grid = size(x)-1;
    L = norm(basis(normalVec+1,:));

    % Correct data so that the sum of all polymer species adds to 1, even
    % in the "wall" area
    for ix = 1:grid(1)+1
        for iy = 1:grid(2)+1
            for iz = 1:grid(3)+1
                pos = [ix,iy,iz];
                d = (pos(normalVec+1)-1) * L / (grid(normalVec+1));
                rho_w = 0.5*(1+tanh(4*(((.5*(T-L))+abs(d-(L/2)))/t)));
                
                if rho_w <= 0.5001
                    R(ix,iy,iz,:) = R(ix,iy,iz,:) ./ (1-rho_w);
                elseif ~ismember(pos(normalVec+1),delete)
                    delete(end+1) = pos(normalVec+1); %#ok<AGROW> 
                end
                
            end
        end
    end
    
    % Delete all data that is "in the wall" (rho_w > 0.5)
    if keep_nan
        if normalVec == 0
            R(delete,:,:,:) = nan;
        elseif normalVec == 1
            R(:,delete,:,:) = nan;
        else % normalVec == 2
            R(:,:,delete,:) = nan;
        end
    else
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
    
end