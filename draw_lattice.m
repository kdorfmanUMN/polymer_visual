% Drawing the Unit Cell Outline:

% This function takes basis, which contains lattice basis vectors, 

function draw_lattice(basis,thick,box_clr)
        
    if nargin == 2 % If box_clr is not provided
        box_clr = [0.5,0.5,0.5];
    elseif nargin == 1 % If box_clr and thick are not provided
        thick = 1;
        box_clr = [0.5,0.5,0.5];
    end

    gca; hold on;

    if size(basis,1) == 3
        % Find coordinates of the corners and store in x_pos, y_pos, and z_pos
        x_pos(1) = 0;
        x_pos(2) = basis(2,1);
        x_pos(3) = basis(1,1) + basis(2,1);
        x_pos(4) = basis(1,1);
        x_pos(5) = basis(3,1);
        x_pos(6) = basis(2,1) + basis(3,1);
        x_pos(7) = basis(1,1) + basis(2,1) + basis(3,1);
        x_pos(8) = basis(1,1) + basis(3,1);
        
        x_start = x_pos;
        for i =9:12
            x_start(i)=x_pos(i-8);
        end
        x_end = zeros(length(x_start),1)';
        
        for i =1:4
            if i == 4
                x_end(i)=x_start(1);
            else
                x_end(i)=x_start(i+1);
            end
        end
        
        for i =5:8
            if i == 8
                x_end(i)=x_start(5);
            else
                x_end(i)=x_start(i+1);
            end
        end
        
        for i =9:12
            x_end(i)=x_start(i-4);
        end
        
        X1 = [x_start;x_end];
        
        y_pos(1) = 0;
        y_pos(2) = basis(2,2);
        y_pos(3) = basis(2,2);
        y_pos(4) = 0;
        y_pos(5) = basis(3,2);
        y_pos(6) = basis(2,2) + basis(3,2);
        y_pos(7) = basis(2,2) + basis(3,2);
        y_pos(8) = basis(3,2);
        
        y_start = y_pos;
        for i =9:12
            y_start(i)=y_pos(i-8);
        end
        y_end = zeros(length(y_start),1)';
        
        for i =1:4
            if i == 4
                y_end(i)=y_start(1);
            else
                y_end(i)=y_start(i+1);
            end
        end
        
        for i =5:8
            if i == 8
                y_end(i)=y_start(5);
            else
                y_end(i)=y_start(i+1);
            end
        end
        
        for i =9:12
            y_end(i)=y_start(i-4);
        end
        
        Y1 = [y_start;y_end];
        
        z_pos(1) = 0;
        z_pos(2) = 0;
        z_pos(3) = 0;
        z_pos(4) = 0;
        z_pos(5) = basis(3,3);
        z_pos(6) = basis(3,3);
        z_pos(7) = basis(3,3);
        z_pos(8) = basis(3,3);
        
        z_start = z_pos;
        for i =9:12
            z_start(i)=z_pos(i-8);
        end
        z_end = zeros(length(z_start),1)';
        
        for i =1:4
            if i == 4
                z_end(i)=z_start(1);
            else
                z_end(i)=z_start(i+1);
            end
        end
        
        for i =5:8
            if i == 8
                z_end(i)=z_start(5);
            else
                z_end(i)=z_start(i+1);
            end
        end
        
        for i =9:12
            z_end(i)=z_start(i-4);
        end
        
        Z1 = [z_start;z_end];
        
        % Plot resulting box
        line(X1,Y1,Z1,'color',box_clr,'LineStyle','-','LineWidth',thick)

    elseif size(basis,1) == 2
        X1 = [0, basis(1,1), basis(1,1)+basis(2,1), basis(2,1)];
        Y1 = [0, 0, basis(2,2), basis(2,2)];
        line(X1,Y1,'color',box_clr,'LineStyle','-','LineWidth',thick)

    else
        error("draw_lattice is only defined for 2D and 3D structures")

    end

end