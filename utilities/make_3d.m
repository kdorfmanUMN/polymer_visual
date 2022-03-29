% Takes 2D rgrid data Ri, as well as coordinate arrays xi and yi, and 
% converts them into a 3D data set [R,x,y,z] with height defined as an
% input variable. We choose arbitrarily to convert our data such that it is
% discretized into 11 z positions from 0 to height (so, 10 gridpoints).

function [R,x,y,z] = make_3d(Ri,xi,yi,height)
    
    R = zeros(size(Ri,1),size(Ri,2),11,size(Ri,3));
    x = zeros(size(Ri(:,:,:,1)));
    y = zeros(size(Ri(:,:,:,1)));
    z = zeros(size(Ri(:,:,:,1)));
    for iz = 1:11
        R(:,:,iz,:) = Ri;
        x(:,:,iz) = xi;
        y(:,:,iz) = yi;
        z(:,:,iz) = (iz-1)*height/10;
    end
    
end