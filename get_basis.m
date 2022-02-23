% Constructs basis vectors in cartesian coordinates based on lattice
% parameters and angles. Generalized for triclinic.

% By convention, align the first axis, with length a, with the x axis.
% The x-y plane is defined to contain lattice vectors a and b. 

% Also, the function has an optional second output, which represents
% the reciprocal lattice vectors corresponding to the real-space
% lattice vectors defined in basis. 

function [basis,kbasis] = get_basis(cell_d,angles)
   
    a = cell_d(1);
    b = cell_d(2);
    c = cell_d(3);
    alpha = angles(1);
    beta = angles(2);
    gamma = angles(3);

    basis(1,1) = a;
    basis(2,1) = b * cos(gamma);
    basis(2,2) = b * sin(gamma);
    basis(3,1) = c * cos(beta);
    basis(3,2) = c * (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma);
    basis(3,3) = c * sqrt(1 - cos(beta)^2 - ((cos(alpha) - cos(beta)*cos(gamma))/sin(gamma))^2);
    basis = round(basis,8);

    kbasis = zeros(size(basis));
    
end