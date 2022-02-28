function polymer_visual(filename)

    close all;

    phase = "C14"; % Leave as "" if you don't want voronoi partition shown on 
                % 2D contour plot. Otherwise, specify phase name
    isovalue = [];

    % Read in the rgrid file
    tic
    [R,x,y,z,dim,lattype,cell_d,angle,n_mnr,grid] = read_rgrid(filename);
    basis = get_basis(cell_d,angle);
    
    comp_disp = 1:n_mnr;
    
    % Other Inputs
    linedraw = true; % draw the isovalue plot
    n_dp = 3; %Number of significant decimal places for color mapping
    
    inputvec = [0 0 1]; %Direction vector for 1-D density plot

    contourvecs = [0 0 0; % Starting corner of contour plot (must be on gridpoint)
                   1 0 0; % Direction of x-axis of contour plot
                   0 1 0];% Direction of y-axis of contour plot

    opacity = ones(2,n_mnr); %Block Opacity
    %opacity = [1,1;0,0.65;1,1];
    
    thick = 1; %Box Thickness Value
    box_clr = [0.5 0.5 0.5]; %Box color

    drawscatter = []; %Block to simulate scattering through
    h_set = 0:3; %Scattering indices
    k_set = 0:1; %Scattering indices
    l_set = 0:1; %Scattering indices
    mono_label = char(1,3);
    for in = 1:n_mnr
        mono_label(in) = char('A'+in-1);
    end

    %weight(1) = 1.2;
    %weight(2) = 1.2;
    %weight(3) = 1.2;

    %map_choice = [2 1 3]; %Map colors
    %mono_disp = [];    %Desired monomers to display
    %comp_disp = [1 3]; %Desired monomers to display in the composite

    max_comps = zeros(1,n_mnr);
    for in = 1:n_mnr
        max_comps(in) = max(max(max(R(:,:,:,in))));
    end
    
    map_store = get_colormaps(n_dp);
    
    % Get isovalues if they are not user-defined:
    if isempty(isovalue)
        isovalue = get_isovalues(R,dim,n_mnr,grid,linedraw,map_store);
    end
    
    % Draw individual density profiles for each monomer species specified
    % in mono_disp:
    disp(ismatrix(R))
    individual_profiles(R,x,y,z,"isovalue",isovalue,...
                        "mono_label",mono_label,"opacity",opacity,...
                        "map",map_store)

    % Drawing the Composite Density Profile
    
    newisovalue = zeros(1,n_mnr);
    cn = zeros(1,n_mnr);
    for in = 1:n_mnr
        newisovalue(in) = isovalue(in) + in - 1;
        %Effective color map range (+5 to buffer)
        cn(in) = 1+ ceil((10^n_dp)*max_comps(in)-((10^n_dp)*isovalue(in))); 
    end
    if ~isempty(comp_disp)
        %
        comp_disp = sort(comp_disp);

        fill = zeros(1,n_mnr-1);

        c = 0;
        for in = comp_disp(1:end-1)
            c = c +1;

            if in~= n_mnr
                fill(in) = (newisovalue(comp_disp(c+1))-newisovalue(in))*(10^n_dp) - cn(in);
            end

        end

        fillmap_store = cell(1,n_mnr-1);
        c = 0;
        for in = comp_disp(1:end-1)
            c = c +1;
            temp_fillmap = zeros(round(fill(in)),3);
            temp_fillmap(:,1) = linspace(map_store{in}(end,1),map_store{in}(1,1),round(fill(in))); %Red
            temp_fillmap(:,2) = linspace(map_store{in}(end,2),map_store{in}(1,2),round(fill(in))); %Green
            temp_fillmap(:,3) = linspace(map_store{in}(end,3),map_store{in}(1,3),round(fill(in))); %Blue
            fillmap_store{in} = temp_fillmap;
        end

        newmap = [];
        for in = comp_disp
            if isovalue(in) < max_comps(in)
                if in == n_mnr
                    newmap = [newmap;map_store(in)];
                else
                    newmap = [newmap;map_store(in);fillmap_store(in)];
                end
            end
        end

        % Simultaneous Visualisation
        if n_mnr > 4
            cdp = figure(); hold on;
            set (cdp, 'Units', 'normalized', 'Position', [0,0,1,1]);
        else
            figure(); hold on;
        end
        title('Composite Density Profile')

        D = R;
        for in = comp_disp
            D(:,:,:,in) = R(:,:,:,in) + in - 1;
            data = (D(:,:,:,in));
            p1 = patch(isosurface(x,y,z,data,newisovalue(in)), ...
                'FaceColor',map_store{in}(1,:),'EdgeColor','none','FaceAlpha',opacity(1,in));
            p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
                'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(2,in));

            if dim == 3
                view(3);                        %Sets the view to 3-D
            else
                view(2);                        %Sets the view to 2-D
            end
            axis equal;                     %Equates the aspect ratio for each axis
            axis vis3d;                     %Freezes aspect ratio (allowing rotation)
            %axis tight;                     %Snaps the axis to the data set

            if strcmp(lattype,'hexagonal') == 1
                size_grid = (grid(1)+1)*(grid(2)+1)*(grid(3)+1);
                coord_set = zeros(size_grid,3);

                % Draw unit cell outlines for the second and third unit cells:
                angle_2 = [pi/2 pi/2 (-2*pi)/3];
                basis_2 = get_basis(cell_d,angle_2);
                draw_lattice(basis_2,thick,box_clr) % Draw second unit cell

                pi3 = -cell_d(1).*cos(pi/3); % Third unit cell in this code block
                pi6 = -cell_d(1).*cos(pi/6);
                line([-1.*cell_d(1);pi3],[0;pi6],[0;0],'color',box_clr,'LineStyle','-',...
                     'LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;-pi6],[0;0],'color',box_clr,'LineStyle','-',...
                     'LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;pi6],[cell_d(3);cell_d(3)],'color',box_clr,...
                     'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;-pi6],[cell_d(3);cell_d(3)],'color',box_clr,...
                     'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);-1.*cell_d(1)],[0;0],[0;cell_d(3)],'color',box_clr,...
                     'LineStyle','-','LineWidth',thick);

                for i =1:2

                    counter = 0;
                    rotangle = 2*pi/3;

                    for iz = 1:grid(3)+1
                        for iy = 1:grid(2)+1
                            for ix = 1:grid(1)+1
                                counter = counter +1;
                                coord_set(counter,1) = x(ix,iy,iz) ;
                                coord_set(counter,2) = y(ix,iy,iz) ;
                                coord_set(counter,3) = z(ix,iy,iz) ;
                            end
                        end
                    end

                    coord_set = coord_set*[cos(rotangle),sin(rotangle),0;-sin(rotangle),cos(rotangle),0;0,0,1];

                    counter = 0;
                    for iz = 1:grid(3)+1
                        for iy = 1:grid(2)+1
                            for ix = 1:grid(1)+1
                                counter = counter +1;
                                x(ix,iy,iz) = coord_set(counter,1) ;
                                y(ix,iy,iz) = coord_set(counter,2) ;
                                z(ix,iy,iz) = coord_set(counter,3) ;
                            end
                        end
                    end

                    data = D(:,:,:,in);
                    p1 = patch(isosurface(x,y,z,data,newisovalue(in)), ...
                        'FaceColor',map_store{in}(1,:),'EdgeColor','none','FaceAlpha',opacity(1,in));
                    p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
                        'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(2,in));
                    set(gcf,'Renderer','zbuffer')
                end
            end

            if isovalue(in) >= max_comps(in)
                % Creating the label for the isovalue
                cell1 = {['\fontsize{14}\Phi' '_' mono_label(in) '='], num2str(round(isovalue(in),2))}; 
                text_disp = {strjoin(cell1)};
                
                % Setting the location and color for the label
                text( x(grid(1)+1,1,round(grid(3)/2)), y(grid(1)+1,1,round(grid(3)/2)), ...
                      z(grid(1)+1,1,round(grid(3)/in)),text_disp,'color',map_store{in}(1,:)) 
            end

        end

        colormap(cell2mat(newmap))

        % cbh= colorbar ;
        % set(cbh,'YTick',0:0.1:3)
        %
        % title(cbh,'\fontsize{14}\Phi')
        % cbh.Location = 'eastoutside';
        % cbh.Color = [1 0 1];

        draw_lattice(basis,thick,box_clr)

        if n_mnr>4
            cbh=colorbar;
            c = 0;
            for in = 1:n_mnr
                c = c+1;
                leftlabel(c) = round(isovalue(in)+in -1,2);
                leftlabel2(c) = round(isovalue(in),2);
                c = c+1;
                leftlabel(c) = round(max_comps(in)+in -1,2);
                leftlabel2(c) = round(max_comps(in),2);
            end

            rightlabel = [0:0.1:1 0.1:0.1:1 0.1:0.1:1];
            h1=axes('position',get(cbh,'position'),'color','none','ylim',...
                    [isovalue(1)-0.01,max_comps(n_mnr)+2.01],'ytick',leftlabel,...
                    'yticklabel',leftlabel2,'xtick',[]);
            set(cbh,'YTick',0:0.1:3,'Yticklabel',rightlabel)
            title(cbh,'\fontsize{14}\Phi')
        else
            if n_mnr < 4
                cb_pos(1,:)= [.68 .11 .04 .815]; %x,y,width,length
                cb_pos(2,:)= [.79 .11 .04 .815]; %x,y,width,length
                cb_pos(3,:)= [.90 .11 .04 .815]; %x,y,width,length
            elseif n_mnr == 4
                cb_pos(1,:)= [.57 .11 .04 .815]; %x,y,width,length
                cb_pos(2,:)= [.68 .11 .04 .815]; %x,y,width,length
                cb_pos(3,:)= [.79 .11 .04 .815]; %x,y,width,length
                cb_pos(4,:)= [.90 .11 .04 .815]; %x,y,width,length
            elseif n_mnr > 4
                cb_pos(1,:)= [.57 .11 .04 .4]; %x,y,width,length
                cb_pos(2,:)= [.68 .11 .04 .4]; %x,y,width,length
                cb_pos(3,:)= [.79 .11 .04 .4]; %x,y,width,length
                cb_pos(4,:)= [.90 .11 .04 .4]; %x,y,width,length
                cb_pos(5,:)= [.57 .11 .04 .4]; %x,y,width,length
                cb_pos(6,:)= [.68 .11 .04 .4]; %x,y,width,length
                cb_pos(7,:)= [.79 .11 .04 .4]; %x,y,width,length
                cb_pos(8,:)= [.90 .11 .04 .4]; %x,y,width,length
            end

            ax(n_mnr+1) = gca; hold on;

            clear cblabel
            for in = 1:n_mnr

                ax(in) = axes;
                colormap(ax(in),cell2mat(map_store(in)))
                ax(in)= gca;
                ax(in).Visible = 'off';
                ax(in).XTick = [];
                ax(in).YTick = [];
                ax(in).ZTick = [];

                if isovalue(in) < max_comps(in)

                    cblabelstart = isovalue(in);
                    cblabelend = max_comps(in);
    %                 if cblabelend - cblabelstart > 0.1
    %                     cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),2);
    %                 else
                    cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),3);
    %                 end
                    l_lngth = linspace(0,1,10);
                    cb(in) = colorbar(ax(in),'Position',cb_pos(in,:));
                    set(cb(in),'ytick',l_lngth,'Yticklabel',cblabel(in,:))
                    %set(cb(in),'ylim',[isovalue(in)-0.01,max_comps(in)+0.01])
                    title1 = strcat('\fontsize{13}\Phi','_',mono_label(in));
                    title(cb(in),title1)
                end
            end

            linkprop(ax,{'view'});
            if n_mnr < 4
                set(ax,'Position',[.05 .11 .55 .815]); %x,y,width,length
            else
                set(ax,'Position',[.05 .11 .44 .815]); %x,y,width,length
            end

            if dim == 3
                view(3);                        %Sets the view to 3-D
            else
                view(2);                        %Sets the view to 2-D
            end
            ax(n_mnr+1) = gca;
            axis vis3d; % Freezes aspect ratio (allowing rotation)
        end
    end

    hold off

    % I(q) plots

    if ~isempty(drawscatter)
        % Scattering

        % Pre-setting the matrix

        count = 0;

        for i = 1:3
            b(i) = (2*pi)/cell_d(i); %basis vectors
        end


        F_store = zeros(length(h_set)*length(k_set)*length(l_set),5);

        for h = h_set
            for k = k_set
                for l = l_set
                    count = count + 1;
                    F_store(count,1) = h;
                    F_store(count,2) = k;
                    F_store(count,3) = l;
                end
            end
        end

        % Calculating the structure factor
        x_index_label = cell(length(F_store),1);
        x_label = cell(length(F_store),1);

        I = zeros(1,length(F_store));
        q = zeros(1,length(F_store));
        x_index = zeros(1,length(F_store));

        for i_f = 1:length(F_store)
            h = F_store(i_f,1);
            k = F_store(i_f,2);
            l = F_store(i_f,3);
            F_sum = 0;
            for in = drawscatter
                x_s = zeros(grid);
                y_s = zeros(grid);
                z_s = zeros(grid);

                for iz=1:grid(3)
                    for iy=1:grid(2)
                        for ix=1:grid(1)
                            x_s(ix,iy,iz) = x(ix,iy,iz)/cell_d(1);
                            y_s(ix,iy,iz) = y(ix,iy,iz)/cell_d(2);
                            z_s(ix,iy,iz) = z(ix,iy,iz)/cell_d(3);
                            F_sum = F_sum + R(ix,iy,iz,in)*exp(2*1i*pi*((h*x_s(ix,iy,iz)) + ...
                                    (k*y_s(ix,iy,iz))+(l*z_s(ix,iy,iz))));

                        end
                    end
                end

            end
            F_store(i_f,4) = F_sum * cell_d(1)*cell_d(2)*cell_d(3);
            F_store(i_f,5) = abs(F_sum)^2;
            I(i_f) = F_store(i_f,5);
            x_index(i_f) = i_f;
            x_index_label(i_f) = {[ h k l ]};
            x_label{i_f}=mat2str(cell2mat(x_index_label(i_f)));
            % x_index_mag(i_f) = h^2 + k^2 + l^2;
            % x_index_3(i_f) = h + k + l;
            q(i_f) = sqrt((b(1)*h)^2 + (b(2)*k)^2 + (b(3)*l)^2);
        end

        % Plotting scatterplots

        %I(q)

        plotmat = sortrows([q;I]');

        [q1,~,q_ind] = uniquetol(plotmat(:,1));

        plotmat_avg = [q1,accumarray(q_ind,plotmat(:,2),[],@mean)];

        q_sort = plotmat_avg(:,1);
        I_sort = plotmat_avg(:,2);

        % Making a delta function

        t_fcr = 100;
        q2 = linspace(0,max(q_sort),length(q_sort)*t_fcr);
        q_plot = sort([q2 q_sort']);
        I_plot = zeros(length(q_plot),1)+min(I);

        step = 1;
        for i = 1:length(q_plot)
            if q_plot(i)== q_sort(step)

                I_plot(i) = I_sort(step);
                step = step + 1;
                if step > length(q_sort)
                    break
                end
            end
        end

        figure(); hold on;
        semilogy(q_plot,I_plot,'.-')
        xlabel(['q [' char(197) '^-^1]'])
        ylabel('Intensity')

        % vs Miller Indices

        t_fcr = 100;
        x_index_plot = linspace(1,length(x_index)*t_fcr,length(x_index)*t_fcr);
        I_plot2 = zeros(length(x_index_plot),1)+min(I);
        for ii = 1:length(x_index)
            I_plot2((ii*t_fcr)-(t_fcr-1)) = I(ii);
        end
        hold off;

        figure(); hold on;
        semilogy(x_index_plot,I_plot2,'.-')
        ylabel('Intensity')
        xlabel('Miller Indices')

        x_label{length(x_label)+1}=' ';

        ax = gca;
        %ax.YLim = [0 max(I)];
        ax.XTick = linspace(1,length(x_index)*t_fcr,length(x_index)+1);
        ax.XTickLabel = x_label;
        ax.XTickLabelRotation = 45;
        
        hold off
    end

    % The 1-D Density Line

    if ~isempty(inputvec)

        startloc = [0 0 0]; %Starting coordinates (normalized) for 1-D density plot

        figure(); hold on
        userinput = inputvec;

        for i = 1:3
            if inputvec(i) < 0
                startloc(i) = startloc(i) + 1;
            end
        end

        start_coord = startloc .* grid + [1 1 1];
        %  end_coord = userinput .* grid + [1 1 1];
        end_coord = start_coord + userinput .* grid;

        dir_vec = end_coord-start_coord;
        step_length = max(abs(dir_vec));
        clear ix iy iz x_plot

        for il = 1:step_length
            ix(il) = start_coord(1)+ round((il-1)*(dir_vec(1)/step_length));
            iy(il) = start_coord(2)+ round((il-1)*(dir_vec(2)/step_length));
            iz(il) = start_coord(3)+ round((il-1)*(dir_vec(3)/step_length));
            x_plot(il) = (il-1)/(step_length-1);
            for in= 1:n_mnr
                line_plot(in,il) = R(ix(il),iy(il),iz(il),in);

            end
        end

        linemarkers(1) = {'s'};
        linemarkers(2) = {'^'};
        linemarkers(3) = {'o'};
        linemarkers(4) = {'p'};
        linemarkers(5) = {'+'};
        linemarkers(6) = {'*'};
        linemarkers(7) = {'d'};
        
        legend_labels3 = strings(1,n_mnr);
        for in = 1:n_mnr
            plot (x_plot,line_plot(in,:),'color',map_store{in}(1,:),...
                  'marker',cell2mat(linemarkers(mod(in,3)+1)),...
                  'markerfacecolor',map_store{in}(1,:),'markersize',5)
            hold on
            legend_labels3(in) = strcat(mono_label(in),' block');

        end


        title1 = strcat('Density Variation Along ','[',num2str(userinput),']');
        title(title1)
        xlabel('\fontsize{14}r/r\fontsize{14}_m_a_x')
        ylabel('\fontsize{14}\Phi\fontsize{13}(r)')
        %ylim([0 1])

        legend(legend_labels3)
        legend('Location','northeastoutside')

        hold off

    end

    % The 2-D Contour Plot

    if ~isempty(contourvecs)

        startloc = contourvecs(1,:); % Starting coordinates (normalized) for contour plot
        xvec = contourvecs(2,:);     % x-axis vector (normalized) for contour plot
        yvec = contourvecs(3,:);     % y-axis vector (normalized) for contour plot

        % Find angle between xvec and yvec by making two new vectors, xvec_orth
        % and yvec_orth that are defined by orthogonal basis vectors
        if or(strcmp(lattype,'tetragonal') == 1, or(strcmp(lattype,'orthorhombic') == 1,...
              strcmp(lattype,'cubic') == 1))
            xvec_orth = xvec;
            yvec_orth = yvec;
        else
            xvec_orth = [xvec(1)+(xvec(2)*cos(angle(3)))+(xvec(3)*cos(angle(2))),...
                         xvec(2)*sin(angle(3)), xvec(3)*sin(angle(2))];
            yvec_orth = [yvec(1)+(yvec(2)*cos(angle(3)))+(yvec(3)*cos(angle(2))),...
                         yvec(2)*sin(angle(3)), yvec(3)*sin(angle(2))];
        end
        cos_theta = dot(xvec_orth,yvec_orth)/(norm(xvec_orth)*norm(yvec_orth));
        theta = acos(cos_theta);

        % Convert vectors to grid coordinates
        start_coord = startloc .* grid + [1 1 1];
        if ~isequal(start_coord,round(start_coord,0))
            disp('Starting coordinates:')
            disp(start_coord)
            error('Starting coordinates are not on a gridpoint')
        end

        xvec_coord = xvec .* grid;
        yvec_coord = yvec .* grid;
        x_step_length = max(abs(xvec_coord));
        y_step_length = max(abs(yvec_coord));
        x_init = linspace(0,norm(xvec),x_step_length);
        y_init = linspace(0,norm(yvec),y_step_length);

        % Find gridpoints & their data vals in the plane of the contour plot
        contour_data = zeros(y_step_length,x_step_length,n_mnr);
        xcontour = zeros(y_step_length,x_step_length);
        ycontour = zeros(y_step_length,x_step_length);
        for xval = 1:x_step_length
            for yval = 1:y_step_length

                x_coord = start_coord(1) + ...
                          (xval-1)*(xvec_coord(1)/x_step_length) + ...
                          (yval-1)*(yvec_coord(1)/y_step_length);
                x_coord = fit_in_cell(x_coord,grid(1));

                y_coord = start_coord(2) + ...
                          (xval-1)*(xvec_coord(2)/x_step_length) + ...
                          (yval-1)*(yvec_coord(2)/y_step_length);
                y_coord = fit_in_cell(y_coord,grid(2));

                z_coord = start_coord(3) + ...
                          (xval-1)*(xvec_coord(3)/x_step_length) + ...
                          (yval-1)*(yvec_coord(3)/y_step_length);
                z_coord = fit_in_cell(z_coord,grid(3));

                %disp([x_coord,y_coord,z_coord])
                if isequal(round([x_coord y_coord z_coord],2),round([x_coord y_coord z_coord],0))
                    contour_data(yval,xval,1:n_mnr) = R(x_coord,y_coord,z_coord,1:n_mnr);
                else
                    contour_data(yval,xval,1:n_mnr) = nan;
                end

                xcontour(yval,xval) = x_init(xval) + (y_init(yval).*cos(theta));
                ycontour(yval,xval) = y_init(yval).*sin(theta);
            end
        end

        % Delete rows/columns if the whole row/column is NaN, because this
        % affects the interpolation of contourf:
        xvals_to_delete = [];
        yvals_to_delete = [];
        for xval = 1:x_step_length
            if isequal(isnan(contour_data(:,xval,1)),ones(size(contour_data(:,xval,1))))
                xvals_to_delete(end+1) = xval;
            end
        end

        for yval = 1:y_step_length
            if isequal(isnan(contour_data(yval,:,1)),ones(size(contour_data(yval,:,1))))
                yvals_to_delete(end+1) = yval;
            end
        end

        contour_data(:,xvals_to_delete,:) = [];
        contour_data(yvals_to_delete,:,:) = [];
        xcontour(:,xvals_to_delete,:) = [];
        xcontour(yvals_to_delete,:,:) = [];
        ycontour(:,xvals_to_delete,:) = [];
        ycontour(yvals_to_delete,:,:) = [];

        
        % Make contour plot
        cblabel2 = zeros(n_mnr,25); 
        cblabel3 = zeros(n_mnr,6);
        
        figure(); hold on; 
        
        for in = 1:n_mnr
            cblabel2(in,:) = round(linspace(isovalue(in),max_comps(in),25),3);
            cblabel3(in,:) = round(linspace(isovalue(in),max_comps(in),6),3);

            if in == 1
                ax(in) = gca;
            else
                ax(in) = axes;
            end

            contourf(xcontour,ycontour,contour_data(:,:,in),cblabel2(in,:),'linecolor','none')
            colormap(gca,cell2mat(map_store(in)))
            cb(in) = colorbar(gca,'Position',cb_pos(in,:));
            set(cb(in),'ytick',cblabel3(in,:),'Yticklabel',cblabel3(in,:))
            caxis([min(cblabel(in,:)),max(cblabel(in,:))])
            title1 = strcat('\phi','_',mono_label(in));
            title(cb(in),title1,'fontsize',20)

            set(ax(in),'Visible','off','fontsize',16,'outerposition',[0 0 .68 1],'fontweight','bold')
            pbaspect(abs( [ norm(xvec_orth.*cell_d.*(1+abs(cos(theta)))), ...
                            norm(yvec_orth.*cell_d.*sin(theta)), 1]))
        end
        set(gcf,'position',[100 100 1200 600]) % 1625 is good on long axis for hexagonal plots
        
        if phase ~= "" % If phase is specified, draw voronoi partition
            
            % Find voronoi cells and boundaries for this crystal
            % structure (all in terms of miller indices!)
            [v,c] = get_voronoi(cell_d,angle,phase);
            n_voro = length(c);
            
            % find lines of voronoi cells intersecting with the plotting plane
            hold on
            for i_vor = 1:n_voro
                ptcloud = v(c{i_vor},:); % corners of voronoi cell

                if any(ptcloud==Inf)
                    continue
                end
                P = get_cross_section(ptcloud,startloc,cross(xvec,yvec));
                if isempty(P)
                    continue
                end
                % find linear combination of plane vectors required to get to these
                % points
                rescaled = ([xvec; yvec]'\P')';
                k = convhull(rescaled);
                plot(rescaled(k,1),rescaled(k,2),'k-','linewidth',2)
            end
            
        end
        xlim([0,1])
        ylim([0,1])

        hold off
        
    end
    toc
end