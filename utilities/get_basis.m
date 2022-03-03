% Constructs basis vectors in Cartesian coordinates based on lattice
% parameters and angles. Generalized for triclinic.

% NOTE: angles for lattice parameters must be in radians, as in PSCF.

% By convention, the first axis, with length a, is aligned with the x axis.
% The x-y plane is defined to contain lattice vectors a and b. 

% The function also has an optional second output, which represents
% the reciprocal lattice vectors corresponding to the real-space
% lattice vectors defined in basis. We use the convention that the dot
% product between a lattice basis vector and its corresponding reciprocal
% lattice basis vector equals 1, because this is the convention used by
% MATLAB's built-in Fast Fourier Transform commands. Rows 1, 2, and 3 of 
% kbasis are the x, y, and z reciprocal lattice vectors, respectively. By
% definition, the dot product between lattice basis vector a_i and
% reciprocal lattice basis vector b_j is zero unless i=j.

function [basis,kbasis] = get_basis(cell_d,angles)
    
    if length(cell_d) == 3 && length(angles) == 3 % 3D
        
        % Pull individual unit cell params from input
        a = cell_d(1);
        b = cell_d(2);
        c = cell_d(3);
        alpha = angles(1);
        beta = angles(2);
        gamma = angles(3);
        
        % Get lattice basis vectors
        basis = zeros(3); 
        basis(1,1) = a;
        basis(2,1) = b * cos(gamma);
        basis(2,2) = b * sin(gamma);
        basis(3,1) = c * cos(beta);
        basis(3,2) = c * (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma);
        basis(3,3) = c * sqrt(1 - cos(beta)^2 - ...
                     ((cos(alpha) - cos(beta)*cos(gamma))/sin(gamma))^2);
        basis = round(basis,10);

        % Get reciprocal lattice basis vectors
        V = abs(dot(basis(1,:),cross(basis(2,:),basis(3,:))));
        kbasis = zeros(3); 
        kbasis(1,:) = cross(basis(2,:),basis(3,:)) / V;
        kbasis(2,:) = cross(basis(3,:),basis(1,:)) / V;
        kbasis(3,:) = cross(basis(1,:),basis(2,:)) / V;
        kbasis = round(kbasis,10);
        
    elseif length(cell_d) == 2 && isscalar(angles) % 2D
        
        % Get lattice basis vectors
        basis = zeros(2);
        basis(1,1) = cell_d(1);
        basis(2,1) = cell_d(2) * cos(angles);
        basis(2,2) = cell_d(2) * sin(angles);
        
        % Get reciprocal lattice basis vectors
        A = basis(1,1) * basis(2,2);
        kbasis = zeros(2); 
        kbasis(1,:) = ([0 1; -1 0] * basis(2,:)' / A)';
        kbasis(2,:) = -1 * ([0 1; -1 0] * basis(1,:)' / A)';
        kbasis = round(kbasis,10);
        
    else % bad input
        
        msg=strcat('cell_d and angles must correspond to a 2D or 3D',...
                   ' unit cell.\nfor 3D, both arrays must have length',...
                   ' 3.\nfor 2D, cell_d must have length 2, and angles',...
                   ' must be %s.');
        error(msg,'scalar')
        
    end
    
end