% This function outputs the 2D aspect ratio of a 3D plot. To do this, it
% accepts az and el, the azimuth and elevation angles that define the
% vector normal to the plane of the figure (the viewing angle). It also
% accepts pb, which should be an array containing the plotbox ratio. Using
% geometry, the function determines and returns the aspect ratio.

% RECOMMENDED USE:
% [az,el] = view;                             % get angles
% pb = pbaspect;                              % get plotbox aspect ratio
% aspect_ratio = get_aspect_ratio(az,el,pb);  % get aspect ratio

% If you are getting anomalous behavior, you may need to add a drawnow()
% preceding the above commands, in order to update your figure to its
% most up-to-date ratios/views.

function aspect_ratio = get_aspect_ratio(az,el,pb)

    if el == 90 % We are looking straight down the z-axis
        
        aspect_ratio = pb(1) / pb(2);
        
    else % We have a 3d view
        
        % view_vec is a vector in data coordinates that points directly at
        % the viewer
        view_vec = [cosd(az-90) sind(az-90) tand(el)];
        view_vec = view_vec / norm(view_vec); % Make it have unit length
        
        % vert_vec is a vector in data coordinates that points in the
        % vertical direction in the plane of the plot (i.e. your screen)
        vert_vec = [0 0 1] - (dot([0 0 1],view_vec) * view_vec);
        vert_vec = vert_vec / norm(vert_vec); % Make it have unit length
        
        % horz_vec is a vector in data coordinates that points in the
        % horizontal direction in the plane of the plot (i.e. your screen)
        horz_vec = cross(vert_vec,view_vec);
        horz_vec = horz_vec / norm(horz_vec); % Make it have unit length
        
        % Get aspect ratio
        width = dot([pb(1) -pb(2) 0], horz_vec); % Reduced plot width
        height = dot(pb,vert_vec); % Reduced plot height
        aspect_ratio = width / height;
        
    end

    aspect_ratio = abs(aspect_ratio);
    
end
