% Drawing a box:

% This function will draw a 3d parallelepiped (hereafter referred to as 
% the "box" for simplicity, though the edges need not be orthogonal) on 
% the current axis. By default, the box edges are drawn in gray with a
% LineWidth of 1, and the box faces are transparent.

% The function requires an input matrix called "basis", the rows of which 
% are the box's basis vectors (the vectors defining the direction of the
% three non-parallel edges of the box). An optional second input called 
% "origin" defines the point in space where the three edges defined in
% "basis" originate (default value is [0 0 0], the actual origin of the
% figure). So, the box will consist of the edges (origin+basis(i,:)) for
% i=1:3, plus the other nine edges needed to complete the parallelepiped.

% Optionally, users can provide Name-Value pair inputs to change the
% appearance of the box, or give the box non-transparent faces. The
% Name-Value pair inputs available and their default values are the 
% following (which are passed directly into the "patch" function, so 
% users can find further information about these inputs in the 
% documentation for "patch"): 
%
%   - LineWidth = 1            (width of lines on box edges)
%   - LineStyle = '-'          (style of lines on box edges)
%   - EdgeColor = [.5,.5,.5]   (color of lines on box edges, default=gray)
%   - EdgeAlpha = 1            (opacity of lines on box edges, from 0-1)
%   - FaceColor = 'none'       (color of box faces)
%   - FaceAlpha = 1            (opacity of box faces, from 0-1)
%   - Marker = 'none'          (marker type for box vertices)
%   - MarkerSize = 6           (size of markers on box vertices)
%   - MarkerEdgeColor = 'auto' (color of edges of markers on box vertices)
%   - MarkerFaceColor = 'none' (color of faces of markers on box vertices)
%
% Users can pass anything into these Name-Value pair arguments as long as
% they can be read by "patch". For example, "EdgeColor" could be a RGB
% triplet like [.5,.5,.5], or it could be a string that indicates a certain
% color such as 'y' for yellow. 

function draw_box(basis,origin,options)
    
    arguments

        % basis must be a 3x3 array where each row is a basis vector.
        basis (3,3) {mustBeNumeric} 

        % origin is a vector indicating the location of the bottom corner
        % of the box to be drawn (i.e., the corner of the box at which the
        % three edges corresponding to the lattice basis vectors
        % originate). This is a 2-element vector in 2d or a 3-element
        % vector in 3d. Default placement of origin is [0 0 0] (which is
        % read by the program as merely [0 0] if basis is 2x2).
        origin (1,3) {mustBeNumeric} = [0 0 0]
        
        % Optional Name-Value pair inputs. See above for details
        options.LineWidth = 1;
        options.LineStyle = '-';
        options.EdgeColor = [0.5,0.5,0.5];
        options.EdgeAlpha = 1;
        options.FaceColor = 'none';
        options.FaceAlpha = 1;
        options.Marker = 'none';
        options.MarkerSize = 6;
        options.MarkerEdgeColor = 'auto';
        options.MarkerFaceColor = 'none';

    end

    % Get current axis and set "hold on" so we can plot over top of any
    % existing graphics objects
    gca; hold on;
    
    % Extract the basis vectors from the input matrix "basis"
    a = basis(1,:);
    b = basis(2,:);
    c = basis(3,:);
    
    % Define the locations of the vertices
    v = [0 0 0; a; a+b; b; c; a+c; a+b+c; b+c]; % 8 vertices

    % Shift all vertices by the optional input "origin"
    for i=1:8
        v(i,:) = v(i,:) + origin;
    end

    % Define the faces based on which vertices they contain
    f = [1 2 3 4 1; 1 5 6 2 1; 1 5 8 4 1; ...
         2 3 7 6 2; 3 4 8 7 3; 5 6 7 8 5]; % 6 faces

    % Create the box using the "patch" function
    patch('Faces', f, 'Vertices', v, ...
          'LineWidth', options.LineWidth, ...
          'LineStyle', options.LineStyle, ...
          'EdgeColor', options.EdgeColor, ...
          'EdgeAlpha', options.EdgeAlpha, ...
          'FaceColor', options.FaceColor, ...
          'FaceAlpha', options.FaceAlpha, ...
          'Marker', options.Marker, ...
          'MarkerSize', options.MarkerSize, ...
          'MarkerEdgeColor', options.MarkerEdgeColor, ...
          'MarkerFaceColor', options.MarkerFaceColor);

end