function polymer_visual(filename)

    close all;

    phase = "C14"; % Leave as "" if you don't want voronoi partition shown on 
                % 2D contour plot. Otherwise, specify phase name
    isovalue = [];

    % Read in the rgrid file
    tic
    [R,x,y,z,dim,lattype,cell_d,angle,n_mnr,grid] = read_rgrid(filename);
    basis = get_basis(cell_d,angle);
    
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
    box_color = [0.5 0.5 0.5]; %Box color
    
    scatterers = 1;
    units = 'Å';
    
    cb_rows = 1; % # of rows of colorbars in composite density profile
    cb_ticks = 8; % # of ticks on colorbars
    n_digits = 3; % # of digits after the decimal point shown in cb ticklabels

    drawscatter = []; %Block to simulate scattering through
    h_set = 0:3; %Scattering indices
    k_set = 0:1; %Scattering indices
    l_set = 0:1; %Scattering indices
    mono_label = char(1,3);
    for in = 1:n_mnr
        mono_label(in) = char('A'+in-1);
    end

    max_comps = zeros(1,n_mnr);
    for in = 1:n_mnr
        max_comps(in) = max(max(max(R(:,:,:,in))));
    end
    
    map_store = get_colormaps();
    fontsize = 15;
    
    % Get isovalues if they are not user-defined:
    if isempty(isovalue)
        isovalue = get_isovalues(R,dim,n_mnr,grid,linedraw,map_store,...
                                 fontsize);
    end
    
    make_3d = false;
    savefile = "";
    writefile = ""; % to write scattering data
    if strcmp(lattype,"hexagonal")
        hex3=true;
    else
        hex3=false;
    end
    
    % Draw individual density profiles for each monomer species specified
    % in mono_disp:
    individual_profiles(R,x,y,z,"isovalue",isovalue,"map",map_store,...
                        "mono_label",mono_label,"opacity",opacity,...
                        "hex3",hex3,"thick",thick,"box_color",box_color,...
                        "make_3d",make_3d,"cb_ticks",cb_ticks,"fontsize",...
                        fontsize,"savefile",savefile)

    % Draw the Composite Density Profile
    composite_profile(R,x,y,z,"isovalue",isovalue,"map",map_store,...
                      "mono_label",mono_label,"opacity",opacity,...
                      "hex3",hex3,"thick",thick,"box_color",box_color,...
                      "make_3d",make_3d,"cb_ticks",cb_ticks,"fontsize",...
                      fontsize,"cb_rows",cb_rows,"n_digits",n_digits,...
                      "savefile",savefile)

    % Draw the scattering plot
    [f_path,f_name,f_ext] = fileparts(savefile);
    scat_savefile = fullfile(f_path,strcat(f_name,'_scat',f_ext));
    scattering_plot(R,x,y,z,scatterers,'fontsize',fontsize,'savefile',...
                    scat_savefile,'theta_plot',false,'units',units,...
                    'writefile',writefile);
    
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
            cb(in) = colorbar(gca);
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