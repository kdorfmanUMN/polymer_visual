function P = get_cross_section(ptcloud,pt,norm)

    % PTCLOUD is a bunch of points falling on the surface of and
    % potentially inside of a 3D volume
    %
    % the volume is assumed to be described by the convex hull of PTCLOUD
    %
    % pt is a point in the plane that is assumed to intersect the 3D
    % volume. 
    %
    % norm is the normal vector of the plane. norm and -norm give the same
    % result.
    %
    % the plane is checked to intersect the ptcloud, and if it does not, 
    % 0 will be returned and a warning displayed.
    %
    % P is the intersections of the plane with edges of the Delaunay
    % triangulation of PTCLOUD
    
    % make sure pt/norm are correct shape
    pt = reshape(pt,1,3);
    norm = reshape(norm,1,3);
    % make sure ptcloud has correct shape
    if size(ptcloud,2)~=3
        disp('ptcloud is wrong size.')
        return
    end
    
    % find convex hull/triangulation of convex hull
    k = convhulln(ptcloud);
    P = [];
    n_tri = size(k,1); % number of triangles on convex hull
    for i_tri = 1:n_tri
        % for each triangle, it is a special case if any one of the corners is in the plane. 
        % If all corners are in the plane, the triangle is coplanar. 
        % If just two are in the plane, one of the lines is in the plane.
        % If just one is in the plane, then that one point touches the plane.
        % In each of these cases, the only point that needs to be added is the corner that touches
        % If none, check if plane intersects lines of triangle
        % If yes, store intersection points as a row in P
        % If no, go to next triangle
        
        % get triangle coords/vectors
        pt_tri = ptcloud(k(i_tri,:),:);
        v = [pt_tri(3,:) - pt_tri(2,:); % vectors of triangle lines
            pt_tri(3,:) - pt_tri(1,:);
            pt_tri(2,:) - pt_tri(1,:)];
        norm_tri = cross(v(2,:),v(3,:)); 
        
        % check if any points of triangle are in plane.
        cornersIn = false;
        for i_corn = 1:3
            if round(dot(norm, pt_tri(i_corn,:)-pt),7) == 0
                P = cat(1,P,pt_tri(i_corn,:));
                cornersIn = true;
            end
        end
        % if any were in, go to the next triangle
        if cornersIn
            continue
        end
        
        % if no points are in the plane, then check if they intersect
        allidx = [1 2 3];
        for i_line = 1:length(allidx)
            % all indexes except for the point not in the line. The point not in line is indexed i_line
            idxs = allidx(allidx~=i_line); 
            % find intersection of line with plane.
            % first, find the value of the scalar multiple of the line vector
            t = dot( norm, (pt-pt_tri(idxs(1),:)) )/dot( norm, v(i_line,:) );
            % calculate the point of intersection using the identified scalar multiple
            pt_int = pt_tri(idxs(1),:)+t*v(i_line,:); 
            % now, check if this intersection of the line is actually on the component of the line in the triangle. 
            % It could be away from the triangle! 
            isIn = isInTri(pt_int,pt_tri);
            if isIn
                P = cat(1,P, pt_int);
            end
            % if not in, go onto the next line.
        end 
        % move onto the next triangle
    end
    % found all points of intersection of plane with triangles
    % return points as an n_pt x 3 matrix
    % draw convex hull of this to get the actual cross-section shape!
end

function bool = isInTri(pt,pt_tri)
    % determines if a point p lies within a triangle (including boundaries)
    % corners of triangle are the rows of pt_tri
    % make sure input is correctly shaped
    pt = reshape(pt,1,3);
    % first, checks if that point is in the plane of the triangle
    v = [pt_tri(3,:) - pt_tri(2,:); % vectors of triangle lines
        pt_tri(3,:) - pt_tri(1,:);
        pt_tri(2,:) - pt_tri(1,:)];
    norm_tri = cross(v(2,:),v(3,:));
    % if not in plane
    if round(dot( norm_tri, pt-pt_tri(1,:) ), 6) ~= 0
        bool = false;
        return
    end
    
    % if it is in plane, then..
    % get linear combo of triangle vectors to get to point
    b = [pt_tri(1,:)' v(3,:)' v(2,:)']\pt';
    
    % if the linear combo of the vectors is within the triangle...
    b = round(b,10);
    if b(2)+b(3) <= 1 && b(2) >= 0 && b(3) >=0 && b(1)==1
        bool = true;
    else
        bool = false;
    end
    
end