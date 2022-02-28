% Drawing the Unit Cell Outline:

% This function will draw the edges of a 2d or 3d box on the current axis. 

% The function takes basis, the rows of which are the lattice basis 
% vectors, and thick and box_clr, which are specifiers for the line
% thickness and box color of the unit cell outline, respectively.
% thick and box_clr are optional; they default to 1 and gray, respectively.

function draw_lattice(basis,thick,box_clr)
    
    arguments

        % basis must be a 2x2 or 3x3 array where each row is a lattice 
        % basis vector, where the array size depends on the dimensionality
        % of the system being plotted.
        basis {mustBeNumeric} 
        
        % thick is the value for "linethickness" used to draw the outer box
        % of the unit cell. Default is 1
        thick = 1;

        % box_color is the value for "color" used to draw the outer box of
        % the unit cell. Default is gray.
        box_clr = [0.5,0.5,0.5];

    end

    gca; hold on;

    if size(basis,1) == 3
        % Find coordinates of the corners and store in x_pos, y_pos, and z_pos
        a = basis(1,:);
        b = basis(2,:);
        c = basis(3,:);

        line_start = zeros(3,12);
        line_end = zeros(3,12);

        line_start(:,2) = a;
        line_start(:,3) = a+b;
        line_start(:,4) = b;
        line_start(:,6) = a;
        line_start(:,7) = a+b;
        line_start(:,8) = b;
        line_start(:,9) = c;
        line_start(:,10) = c+a;
        line_start(:,11) = c+a+b;
        line_start(:,12) = c+b;

        line_end(:,1) = a;
        line_end(:,2) = a+b;
        line_end(:,3) = b;
        line_end(:,5) = c;
        line_end(:,6) = c+a;
        line_end(:,7) = c+a+b;
        line_end(:,8) = c+b;
        line_end(:,9) = c+a;
        line_end(:,10) = c+a+b;
        line_end(:,11) = c+b;
        line_end(:,12) = c;

        X1 = [line_start(1,:); line_end(1,:)];
        Y1 = [line_start(2,:); line_end(2,:)];
        Z1 = [line_start(3,:); line_end(3,:)];

        % Plot resulting box
        line(X1,Y1,Z1,'color',box_clr,'LineStyle','-','LineWidth',thick)

    elseif size(basis,1) == 2
        X1 = [0, basis(1,1), basis(1,1)+basis(2,1), basis(2,1)];
        Y1 = [0, basis(1,2), basis(1,2)+basis(2,2), basis(2,2)];
        line(X1,Y1,'color',box_clr,'LineStyle','-','LineWidth',thick)

    else
        error("draw_lattice is only defined for 2D and 3D structures")

    end

end