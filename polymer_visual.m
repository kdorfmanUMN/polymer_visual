function polymer_visual(phase,filename)

    close all

    % Read in the rgrid file
    tic
    [R,x,y,z,dim,lattype,cell_d,angle,n_mnr,grid] = read_rgrid(filename);

    % Default Inputs
    mono_disp = zeros(1,n_mnr);
    map_choice= zeros(1,n_mnr);
    comp_disp = zeros(1,n_mnr);
    weight = zeros(1,n_mnr);
    for in = 1:n_mnr
        mono_disp(in) = in;
        map_choice(in)= in;
        comp_disp(in) = in;
        weight(in) = 1;
    end
    
    % Other Inputs

    ncolor = 8; %Number of stored color maps
    n_dp = 3; %Number of significant decimal places for color mapping
    gap_close = 1; %Toggle the white space closing feature (Default = 1)

    inputvec = [0 0 1]; %Direction vector for 1-D density plot
    % C15:
    % contourvecs = [0 0 0; % Starting corner of contour plot (must be on gridpoint)
    %               1 1 0; % Direction of x-axis of contour plot
    %               0 0 1];% Direction of y-axis of contour plot
    % C14:
    % contourvecs = [1/3 -1/3 7/16; % Starting corner of contour plot (must be on gridpoint)
    %                2 0 0; % Direction of x-axis of contour plot
    %                0 2 0];% Direction of y-axis of contour plot
    % sigma:
    contourvecs = [0 0 0; % Starting corner of contour plot (must be on gridpoint)
                  1 0 0; % Direction of x-axis of contour plot
                  0 1 0];% Direction of y-axis of contour plot

    opacity = ones(n_mnr,2); %Block Opacity
    %opacity = [1,1;0,0.65;1,1];
    isovalue = 'auto';
    thick = 1; %Box Thickness Value
    box_clr = [0.5 0.5 0.5]; %Box color

    drawscatter = []; %Block to simulate scattering through
    h_set = 0:3; %Scattering indices
    k_set = 0:1; %Scattering indices
    l_set = 0:1; %Scattering indices

    %weight(1) = 1.2;
    %weight(2) = 1.2;
    %weight(3) = 1.2;

    %map_choice = [2 1 3]; %Map colors
    %mono_disp = [];    %Desired monomers to display
    %comp_disp = [1 3]; %Desired monomers to display in the composite

    % Reading Additional Inputs From A Text File (input.txt)

    % C2 = textread('input.txt', '%s','delimiter', '\n');
    %
    % ic = 0;
    %
    % for ic = 1:length(C2)
    %
    %     if strcmp(strrep(char(C2(ic)),' ', ''),'thickness')==1
    %         thick = str2double(C2{ic+1});          % Reads the box thickness
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'box_clr')==1
    %         box_clr = sscanf(C2{ic+1},'%f');           % Reads the box color
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'opacity')==1
    %         opacity(:,1) = sscanf(C2{ic+1},'%f')'; % Reads the opacity
    %         opacity(:,2) = sscanf(C2{ic+2},'%f')'; % Reads the opacity
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'isovalue')==1
    %         if ischar(strrep(char(C2(ic+1)),' ', ''))
    %         isovalue = strrep(char(C2(ic+1)),' ', '');
    %         else
    %         isovalue = sscanf(C2{ic+1},'%f')';           % Reads the isovalues
    %         end
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'map_choice')==1
    %         map_choice = sscanf(C2{ic+1},'%f')';    % Reads the chosen map colors
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'mono_disp')==1
    %         mono_disp = sscanf(C2{ic+1},'%f')';           % Reads the desired monomers to display
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'comp_disp')==1
    %         comp_disp = sscanf(C2{ic+1},'%f')';          % Reads the desired monomers to display in the composite
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'weight')==1
    %         weight = sscanf(C2{ic+1},'%f')';          % Reads the desired weight of composite monomers
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'scatterdraw')==1
    %         drawscatter = sscanf(C2{ic+1},'%f')';     % For scattering
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'h_set')==1
    %         h_set = sscanf(C2{ic+1},'%f')';     % For scattering
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'k_set')==1
    %         k_set = sscanf(C2{ic+1},'%f')';     % For scattering
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'l_set')==1
    %         l_set = sscanf(C2{ic+1},'%f')';     % For scattering
    %     elseif strcmp(strrep(char(C2(ic)),' ', ''),'input_vec')==1
    %         input_vec = sscanf(C2{ic+1},'%f')';          %For the 1-D Plot
    %     end
    % end

    if strcmp(isovalue,'auto')==1
        linedraw = 1; %Diagnostic Plot
    else
        linedraw=0;
    end
    %linedraw=0; %(Set to 0 if you do not want to view the diagnostic plot)

    % Computing the Isovalue

    for ig = 1:3
        plot_grid(ig) = grid(ig)+1;
    end

    if strcmp(isovalue,'auto')==1

        isovalue = zeros(1,n_mnr);

        %Finding the range and limits
        n_lines = 1;

        polmaxa = zeros(1,n_mnr);
        polmin = zeros(1,n_mnr);
        l_length = zeros(1,n_mnr);
        for in = 1:n_mnr
            polmaxa(in) = max(max(max((R(:,:,:,in)))));
            polmin(in) = min(min(min((R(:,:,:,in)))));
            l_length(in) = polmaxa(in) -polmin(in);
        end


        %Creating the scaled matrix, S
        S = R;
        for in = 1:n_mnr
            for iz=1:grid(3)+1
                for iy=1:grid(2)+1
                    for ix=1:grid(1)+1
                        S (ix,iy,iz,in) = weight(in)*(R(ix,iy,iz,in)-polmin(in))/l_length(in); % Matrix with scaled values
                    end
                end
            end
        end

        % Finding the location of the polmax points
        pmax_loc = zeros(n_lines,3,n_mnr);
        for in =1:n_mnr
            [si,sj,sk] = ind2sub([plot_grid n_mnr],find(S(:,:,:,in) == weight(in),n_lines));
            pmax_loc(:,:,in) = [si,sj,sk];
        end

        % Interleaving
        point_series = zeros(n_lines*n_mnr,3);
        for ir = 1:n_lines
            for in = 1:n_mnr
                point_series(n_mnr*(ir-1)+in,:) = pmax_loc(ir,:,in); %Series of points to be plotted
            end
        end

        % Creating the lines
        clear line line_new x_plot intervalue ix iy iz
        dir_vec = zeros((n_lines*n_mnr)-1,3);
        step_length = zeros((n_lines*n_mnr)-1);

        ix = zeros((n_lines*n_mnr)-1,step_length(ir));
        iy = zeros((n_lines*n_mnr)-1,step_length(ir));
        iz = zeros((n_lines*n_mnr)-1,step_length(ir));
        counter = 0;
        for ir = 1:(n_lines*n_mnr)-1 %Number of lines
            start_coord = point_series(ir,:);
            end_coord = point_series(ir+1,:);
            dir_vec(ir,:) = end_coord-start_coord;
            step_length(ir) = max(abs(dir_vec(ir,:)));

            for il = 1:step_length(ir)
                ix(ir,il) = start_coord(1)+ round((il-1)*(dir_vec(ir,1)/step_length(ir)));
                iy(ir,il) = start_coord(2)+ round((il-1)*(dir_vec(ir,2)/step_length(ir)));
                iz(ir,il) = start_coord(3)+ round((il-1)*(dir_vec(ir,3)/step_length(ir)));
                counter = counter + 1;
                for in= 1:n_mnr
                    line_vals(in,counter) = R(ix(ir,il),iy(ir,il),iz(ir,il),in);
                end
            end
        end
        x_plot = linspace(1,counter+1,counter+1);

        % The final point
        for in= 1:n_mnr
            line_vals(in,counter+1)= R(end_coord(1),end_coord(2),end_coord(3),in);
        end

        % Rescaling the line
        line_new = zeros(n_mnr,length(x_plot));
        for in = 1:n_mnr
            for ip = 1:length(x_plot)
                line_new(in,ip) = weight(in)*(line_vals(in,ip)-polmin(in))/l_length(in);
            end
        end

    % Finding intersections

        clear x_inter inter_point intervalue

        n_intervalue  = 0;
        for in = 1:n_mnr
            n_intervalue  = n_intervalue  + (in-1);
        end

        x_inter_store = cell(1,n_intervalue);
        loop_round = 1;
        k = 0;
        involved = zeros(n_mnr,n_mnr-1);

        while n_mnr-loop_round > 0

            for j = 1:n_mnr-loop_round
                k = k+1;
                involved(loop_round,j)= k;  %Needed for isovalue calculation
                involved(loop_round+j,n_mnr-j)= k;
                % Check if the lines intersect:
                if ~isempty(find(diff(sign(line_new(loop_round+j,:)-line_new(loop_round,:))), 1))
                    x_inter_store{k} = find(diff(sign(line_new(loop_round+j,:)-line_new(loop_round,:))));

                    % For each intersection (there may be multiple),
                    % interpolate to find the y value where it intersects
                    for i = 1:length(cell2mat(x_inter_store(k)))
                        x_inter = cell2mat(x_inter_store(k));
                        p(1,1)= line_new(loop_round,x_inter(i));
                        p(1,2)= line_new(loop_round,1+x_inter(i));
                        p(2,1)= line_new(loop_round+j,x_inter(i));
                        p(2,2)= line_new(loop_round+j,1+ x_inter(i));
                        inter_point(k,i) = (p(1,1)*(p(2,2)-p(2,1))-p(2,1)*(p(1,2)-p(1,1))) / ...
                                           (p(2,2)-p(2,1)-p(1,2)+p(1,1));
                    end
                    ipts = inter_point(k,:);
                    if ~isempty(ipts(ipts > 0.05 & ipts< 0.95))
                        intervalue(k) = max(ipts(ipts < 0.95));
                    else
                        intervalue(k) = max(inter_point(k,:));
                    end
                end
            end
            loop_round = loop_round +1;
        end

        in_mat = zeros(n_mnr,n_mnr-1);
        isovalue_s = zeros(1,n_mnr);
        for k = 1:n_mnr-1
            for in = 1:n_mnr
                in_mat(in,k)= intervalue(involved(in,k));
                if 1 - max(in_mat(in,:)) > 1e-3
                    isovalue_s(in) = max(in_mat(in,:));
                else 
                    isovalue_s(in) = min(in_mat(in,:));
                end

            end
        end

    % Closing any gaps

        if gap_close == 1

            if n_mnr ==3 && range(isovalue_s) > 1e-4

                no_2 = zeros(1,n_mnr);
                for in = 1:n_mnr
                    r_row = in_mat(in,:);
                    no_2(in) = max(r_row(r_row <max(r_row)));
                end

                for in = 1:n_mnr-1
                    if isovalue_s(in)== no_2(in+1)
                        start_mat = find(diff(sign(isovalue_s(in+1)-line_new(in+1,:))));
                        intersect_val = (find(diff(sign(line_new(in,:)-line_new(in+1,:))), 1));
                        start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));
                        py1= line_new(in+1,start_ind);
                        py2= line_new(in+1,1+start_ind);
                        px1= start_ind;
                        px2= start_ind+1;
                        y_int = isovalue_s(in+1);
                        x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);


                        py1= line_new(in,start_ind);
                        py2= line_new(in,1+start_ind);
                        px1= start_ind;
                        px2= start_ind+1;
                        y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                        isovalue_s(in) = y_new;
                    end
                end

                for in = 2:n_mnr
                    if isovalue_s(in)== no_2(in-1)
                        start_mat = find(diff(sign(isovalue_s(in-1)-line_new(in-1,:))));
                        intersect_val = (find(diff(sign(line_new(in,:)-line_new(in-1,:))), 1));
                        start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));

                        py1= line_new(in-1,start_ind);
                        py2= line_new(in-1,1+start_ind);
                        px1= start_ind;

                        px2= start_ind+1;
                        y_int = isovalue_s(in-1);
                        x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);


                        py1= line_new(in,start_ind);
                        py2= line_new(in,1+start_ind);
                        px1= start_ind;
                        px2= start_ind+1;
                        y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                        isovalue_s(in) = y_new;
                    end
                end

                if isovalue_s(n_mnr) == no_2(1)
                    start_mat = find(diff(sign(isovalue_s(1)-line_new(1,:))));
                    intersect_val = (find(diff(sign(line_new(n_mnr,:)-line_new(1,:))), 1));
                    start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));

                    py1= line_new(1,start_ind);
                    py2= line_new(1,1+start_ind);
                    px1= start_ind;

                    px2= start_ind+1;
                    y_int = isovalue_s(1);
                    x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);


                    py1= line_new(n_mnr,start_ind);
                    py2= line_new(n_mnr,1+start_ind);
                    px1= start_ind;
                    px2= start_ind+1;
                    y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                    isovalue_s(n_mnr) = y_new;
                end

                if isovalue_s(1) == no_2(n_mnr)
                    start_mat = find(diff(sign(isovalue_s(n_mnr)-line_new(n_mnr,:))));
                    intersect_val = (find(diff(sign(line_new(1,:)-line_new(n_mnr,:))), 1));
                    start_ind = start_mat(abs(start_mat - intersect_val)==min(abs(start_mat - intersect_val)));

                    py1= line_new(n_mnr,start_ind);
                    py2= line_new(n_mnr,1+start_ind);
                    px1= start_ind;
                    px2= start_ind+1;
                    y_int = isovalue_s(n_mnr);
                    x_int = px1 + ((y_int-py1)*(px2-px1))/(py2-py1);


                    py1= line_new(1,start_ind);
                    py2= line_new(1,1+start_ind);
                    px1= start_ind;
                    px2= start_ind+1;
                    y_new = py1 + ((x_int-px1)*(py2-py1))/(px2-px1);
                    isovalue_s(1) = y_new;
                end
            end

        end

    % Converting the isovalue

        for in=1:n_mnr
            isovalue(in) = (isovalue_s(in)*l_length(in))/weight(in) +polmin(in);
        end

    end

    % Drawing the Color Maps

    polmaxf = zeros(1,n_mnr);
    cn = zeros(1,ncolor);

    for j = 1:ncolor
        cn(j)=(10^n_dp);
    end

    face1 = zeros(size(R,2)*size(R,3),n_mnr);
    face2 = zeros(size(R,1)*size(R,3),n_mnr);
    face3 = zeros(size(R,1)*size(R,2),n_mnr);
    face4 = zeros(size(face1));
    face5 = zeros(size(face2));
    face6 = zeros(size(face3));
    face_data = zeros((size(face1,1)*2)+(size(face2,1)*2)+(size(face3,1)*2),n_mnr);
    newisovalue = zeros(1,n_mnr);
    mono_label = strings(1,n_mnr);
    titles = strings(1,n_mnr);
    for in = 1:n_mnr
        face1(:,in) = reshape(squeeze(R(1,:,:,in)),[],1);
        face2(:,in) = reshape(squeeze(R(:,1,:,in)),[],1);
        face3(:,in) = reshape(squeeze(R(:,:,1,in)),[],1);
        face4(:,in) = reshape(squeeze(R(grid(1)+1,:,:,in)),[],1);
        face5(:,in) = reshape(squeeze(R(:,grid(2)+1,:,in)),[],1);
        face6(:,in) = reshape(squeeze(R(:,:,grid(3)+1,in)),[],1);
        face_data(:,in) = [face1(:,in); face2(:,in); face3(:,in); face4(:,in); face5(:,in); face6(:,in)];
        polmaxf(in) = max(face_data(:,in)); %Max Density value on face for polymer i
        cn(map_choice(in)) = 1+ ceil((10^n_dp)*polmaxf(in)-((10^n_dp)*isovalue(in))); %Effective color map range (+5 to buffer)
        newisovalue(in) = in + isovalue(in) - 1;
        mono_label(in) = char(in+'A'-1);
        titles(in) = strcat(mono_label(in),' block Density Profile');
    end


    % low = low fraction, i.e. light

    color_low = zeros (ncolor,3);
    color_low(1,:) = [0,0.7,0.9]; %blue
    color_low(2,:) = [0.9,0,0]; %red
    color_low(3,:) = [1,1,0]; %yellow
    color_low(4,:) = [0,0.9,0.2]; %green
    color_low(5,:) = [0.5,0,1]; %purple
    color_low(6,:) = [1,0,1]; %pink
    color_low(7,:) = [1,0.5,0]; %orange
    color_low(8,:) = [0.75,0.75,0.75]; %grey

    % high = high fraction, i.e. dark

    color_high = zeros (ncolor,3);
    color_high(1,:) = [0,0,0.4]; %dark blue
    color_high(2,:) = [0.4,0,0]; %dark red
    color_high(3,:) = [0.4,0.4,0]; %dark yellow
    color_high(4,:) = [0,0.4,0]; %dark green
    color_high(5,:) = [0.15,0,0.30]; %dark purple
    color_high(6,:) = [0.25,0,0.25]; %dark pink
    color_high(7,:) = [0.3216, 0.1882, 0.1882]; % dark "orange" (brown)
    color_high(8,:) = [0.3,0.3,0.3]; %dark grey


    colorpad(:,:,1) = color_low;
    colorpad(:,:,2) = color_high;
    outcolor = colorpad(:,:,1);

    map_store = cell(1,ncolor);

    for in = 1:ncolor
        temp_map = zeros(cn(in),3);
        temp_map(:,1) = linspace(colorpad(in,1,1),colorpad(in,1,2),cn(in)); %Red
        temp_map(:,2) = linspace(colorpad(in,2,1),colorpad(in,2,2),cn(in)); %Green
        temp_map(:,3) = linspace(colorpad(in,3,1),colorpad(in,3,2),cn(in)); %Blue
        map_store{in}=temp_map;
    end

    % Drawing the Individual Profiles

    for in = mono_disp

        figure(in)
        title(titles(in))
        data = R(:,:,:,in);
        patch(isosurface(x,y,z,data,isovalue(in)), ...
              'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none',...
              'FaceAlpha',opacity(in,1));
        patch(isocaps(x,y,z,data,isovalue(in)), ...
              'FaceColor','interp','EdgeColor','none',...
              'FaceAlpha',opacity(in,2));
        colormap(cell2mat(map_store(map_choice(in))))

        if isovalue(in) < max(face_data(:,in)) %Deciding which colorbars are neccesary

            cblabelstart = isovalue(in);
            cblabelend = polmaxa(in);
            if cblabelend - cblabelstart > 0.1
                cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),2);
            else
                cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),3);
            end
            l_lngth = linspace(cblabelstart,cblabelend,10);
            cbh = colorbar;
            set(cbh,'ylim',[cblabelstart cblabelend],'ytick',...
                l_lngth,'Yticklabel',cblabel(in,:))

            title1 = strcat('\fontsize{13}\Phi','_',mono_label(in));
            title(cbh,title1)


        else
            text_disp = strcat('\fontsize{14}\Phi','_',mono_label(in),'=',...
                               num2str(round(isovalue(in),2))); % Creating the label for the isovalue
            %[s1,s2,s3] = ind2sub([plot_grid n_mnr],find(abs((R(:,:,:,in) - isovalue(in)))<0.001,1));
            text(x(grid(1)+1,1,round(grid(3)/2)), ...
                 y(grid(1)+1,1,round(grid(3)/2)), ...
                 z(grid(1)+1,1,round(grid(3)/in)),...
                 text_disp,'color',outcolor(in,:)) % Setting the location for the label
        end

        if strcmp(type,'hexagonal') == 1

            % Draw unit cell outlines for the second and third unit cells:
            angle_2 = [pi/2 pi/2 (-2*pi)/3];
            draw_lattice(cell_d,angle_2,thick,box_clr) % Draw second unit cell

            pi3 = -cell_d(1).*cos(pi/3); % Third unit cell in this code block
            pi6 = -cell_d(1).*cos(pi/6);

            line([-1.*cell_d(1);pi3],[0;pi6],[0;0],'color',box_clr,'LineStyle','-','LineWidth',thick);
            line([-1.*cell_d(1);pi3],[0;-pi6],[0;0],'color',box_clr,'LineStyle','-','LineWidth',thick);
            line([-1.*cell_d(1);pi3],[0;pi6],[cell_d(3);cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);
            line([-1.*cell_d(1);pi3],[0;-pi6],[cell_d(3);cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);
            line([-1.*cell_d(1);-1.*cell_d(1)],[0;0],[0;cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);

            for i =1:2
                size_grid = (grid(1)+1)*(grid(2)+1)*(grid(3)+1);
                coord_set = zeros(size_grid,3);
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

                figure(in)
                data = R(:,:,:,in);
                p1 = patch(isosurface(x,y,z,data,isovalue(in)), ...
                    'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
                p2 = patch(isocaps(x,y,z,data,isovalue(in)), ...
                    'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));

            end
        end

        if dim == 3
            view(3);                        %Sets 3-D view
        elseif dim == 2
            view(2);                        %Sets 2-D view
        elseif dim == 1
            view(2);                        %Sets 2-D view
        end

        axis equal;                     %Equates the aspect ratio for each axis
        axis vis3d;                     %Freezes aspect ratio (allowing rotation)
        %axis tight;                     %Snaps the axis to the data set

        draw_lattice(cell_d,angle,thick,box_clr)

        set(gcf,'Renderer','zbuffer');

        rotate3d

        % Lighting options
        %   set(p2,'AmbientStrength',.6);
        %   set(p1,'AmbientStrength',.5);

        %   isonormals(data,p1);
        %   lightangle(45,30);
        %   lighting phong

    end

    % Drawing the Composite Density Profile

    if ~isempty(comp_disp)
        comp_disp = sort(comp_disp);

        fill = zeros(1,n_mnr-1);

        c = 0;
        for in = comp_disp(1:end-1)
            c = c +1;

            if in~= n_mnr
                fill(in) = (newisovalue(comp_disp(c+1))-newisovalue(in))*(10^n_dp) - cn(map_choice(in));
            end

        end

        fillmap_store = cell(1,n_mnr-1);
        c = 0;
        for in = comp_disp(1:end-1)
            c = c +1;
            temp_fillmap = zeros(round(fill(in)),3);
            temp_fillmap(:,1) = linspace(colorpad(map_choice(comp_disp(c)),1,2),colorpad(map_choice(comp_disp(c+1)),1,1),round(fill(in))); %Red
            temp_fillmap(:,2) = linspace(colorpad(map_choice(comp_disp(c)),2,2),colorpad(map_choice(comp_disp(c+1)),2,1),round(fill(in))); %Green
            temp_fillmap(:,3) = linspace(colorpad(map_choice(comp_disp(c)),3,2),colorpad(map_choice(comp_disp(c+1)),3,1),round(fill(in))); %Blue
            fillmap_store{in} = temp_fillmap;
        end

        % for in = comp_disp(1:end-1)
        %    temp_fillmap = zeros(round(fill(in)),3);
        %     for ifl =1:fill(in)
        %         temp_fillmap(ifl,:) = temp_map(ifl,:);
        %     end
        %     fillmap_store{in} = temp_fillmap;
        % end


        newmap = [];
        for in = comp_disp
            if isovalue(in) < max(face_data(:,in))
                if in == n_mnr
                    newmap = [newmap;map_store(map_choice(in))];
                else
                    newmap = [newmap;map_store(map_choice(in));fillmap_store(in)];
                end
            end
        end

        % Simultaneous Visualisation
        D = R;
        figure(n_mnr+1)
        title('Composite Density Profile')
        if n_mnr > 4
            cdp = figure(n_mnr+1);
            set (cdp, 'Units', 'normalized', 'Position', [0,0,1,1]);
        end
        for in = comp_disp
            D(:,:,:,in) = R(:,:,:,in) +in -1;
            data = (D(:,:,:,in));
            p1 = patch(isosurface(x,y,z,data,newisovalue(in)), ...
                'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
            p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
                'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));

            if dim == 3
                view(3);                        %Sets the view to 3-D
            else
                view(2);                        %Sets the view to 2-D
            end
            axis equal;                     %Equates the aspect ratio for each axis
            axis vis3d;                     %Freezes aspect ratio (allowing rotation)
            %axis tight;                     %Snaps the axis to the data set

            if strcmp(type,'hexagonal') == 1
                size_grid = (grid(1)+1)*(grid(2)+1)*(grid(3)+1);
                coord_set = zeros(size_grid,3);

                % Draw unit cell outlines for the second and third unit cells:
                angle_2 = [pi/2 pi/2 (-2*pi)/3];
                draw_lattice(cell_d,angle_2,thick,box_clr) % Draw second unit cell

                pi3 = -cell_d(1).*cos(pi/3); % Third unit cell in this code block
                pi6 = -cell_d(1).*cos(pi/6);
                line([-1.*cell_d(1);pi3],[0;pi6],[0;0],'color',box_clr,'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;-pi6],[0;0],'color',box_clr,'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;pi6],[cell_d(3);cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);pi3],[0;-pi6],[cell_d(3);cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);
                line([-1.*cell_d(1);-1.*cell_d(1)],[0;0],[0;cell_d(3)],'color',box_clr,'LineStyle','-','LineWidth',thick);

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
                        'FaceColor',outcolor(map_choice(in),:),'EdgeColor','none','FaceAlpha',opacity(in,1));
                    p2 = patch(isocaps(x,y,z,data,newisovalue(in)), ...
                        'FaceColor','interp','EdgeColor','none','FaceAlpha',opacity(in,2));
                    set(gcf,'Renderer','zbuffer')
                end
            end

            if isovalue(in) >= max(face_data(:,in))
                cell1 = {['\fontsize{14}\Phi' '_' mono_label(in) '='], num2str(round(isovalue(in),2))}; % Creating the label for the isovalue
                text_disp = {strjoin(cell1)};
                text( x(grid(1)+1,1,round(grid(3)/2)), y(grid(1)+1,1,round(grid(3)/2)), z(grid(1)+1,1,round(grid(3)/in)),[text_disp],'color',outcolor(in,:)) % Setting the location and color for the label
            end

        end


        colormap(cell2mat(newmap))

        % cbh= colorbar ;
        % set(cbh,'YTick',0:0.1:3)
        %
        % title(cbh,'\fontsize{14}\Phi')
        % cbh.Location = 'eastoutside';
        % cbh.Color = [1 0 1];

        draw_lattice(cell_d,angle,thick,box_clr)

        if n_mnr>4
            figure(n_mnr+1)
            cbh=colorbar;
            c = 0;
            for in = 1:n_mnr
                c = c+1;
                leftlabel(c) = round(isovalue(in)+in -1,2);
                leftlabel2(c) = round(isovalue(in),2);
                c = c+1;
                leftlabel(c) = round(polmaxa(in)+in -1,2);
                leftlabel2(c) = round(polmaxa(in),2);
            end

            rightlabel = [0:0.1:1 0.1:0.1:1 0.1:0.1:1];
            h1=axes('position',get(cbh,'position'),'color','none','ylim',...
                    [isovalue(1)-0.01,polmaxa(n_mnr)+2.01],'ytick',leftlabel,...
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

            figure(n_mnr+1)
            ax(n_mnr+1) = gca;

            clear cblabel
            for in = 1:n_mnr
                c=0;

                figure(n_mnr+1)

                ax(in) = axes;
                colormap(ax(in),cell2mat(map_store(map_choice(in))))
                ax(in)= gca;
                ax(in).Visible = 'off';
                ax(in).XTick = [];
                ax(in).YTick = [];
                ax(in).ZTick = [];

                if isovalue(in) < max(face_data(:,in))

                    cblabelstart = isovalue(in);
                    cblabelend = polmaxa(in);
    %                 if cblabelend - cblabelstart > 0.1
    %                     cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),2);
    %                 else
                    cblabel(in,:) = round(linspace(cblabelstart,cblabelend,10),3);
    %                 end
                    l_lngth = linspace(0,1,10);
                    cb(in) = colorbar(ax(in),'Position',cb_pos(in,:));
                    set(cb(in),'ytick',l_lngth,'Yticklabel',cblabel(in,:))
                    %set(cb(in),'ylim',[isovalue(in)-0.01,polmaxa(in)+0.01])
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

    % Plotting the rescaled volume fractions and isovalues

    if linedraw == 1
        figure(n_mnr+ 2)

        linestyles(1) = {'-.'};
        linestyles(2) = {':'};
        linestyles(3) = {'--'};
        for in = 1:n_mnr
            plot (x_plot,line_new(in,:),'color',outcolor(map_choice(in),:))
            hold on
        end

        clear yy
        yy = ones(n_mnr,length(x_plot));
        for in = 1:n_mnr
            yy(in,:) = yy(in,:)*isovalue_s(in);
        end

        for in = 1:n_mnr
            figure (n_mnr+2)
            hold on
            plot(x_plot,yy(in,:),'color',outcolor(map_choice(in),:),...
                 'LineWidth',1+n_mnr-in,'LineStyle',...
                 cell2mat(linestyles(mod(in,3)+1)))
        end

        xlabel('Line Index')
        ylabel('Rescaled Volume Fraction')

        legend_labels1 = strings(1,n_mnr);
        legend_labels2 = strings(1,n_mnr);
        for in=1:n_mnr
            mono_label(in) = char(in+'A'-1);
            legend_labels1(in) = strcat(mono_label(in),' block Rescaled Density');
            legend_labels2(in) = strcat(mono_label(in),' block Rescaled Isovalue');
        end


        legend([legend_labels1 legend_labels2])
        legend('Location','best')

    end

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

        figure(n_mnr+3)
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

        figure(n_mnr+4)
        semilogy(x_index_plot,I_plot2,'.-')
        ylabel('Intensity')
        xlabel('Miller Indices')

        x_label{length(x_label)+1}=' ';

        ax = gca;
        %ax.YLim = [0 max(I)];
        ax.XTick = linspace(1,length(x_index)*t_fcr,length(x_index)+1);
        ax.XTickLabel = x_label;
        ax.XTickLabelRotation = 45;

    end

    % The 1-D Density Line

    if ~isempty(inputvec)

        startloc = [0 0 0]; %Starting coordinates (normalized) for 1-D density plot

        figure(n_mnr+5)
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

        for in = 1:n_mnr
            plot (x_plot,line_plot(in,:),'color',outcolor(map_choice(in),:),...
                  'marker',cell2mat(linemarkers(mod(in,3)+1)),...
                  'markerfacecolor',outcolor(map_choice(in),:),'markersize',5)
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



    end

    % The 2-D Contour Plot

    if ~isempty(contourvecs)

        startloc = contourvecs(1,:); % Starting coordinates (normalized) for contour plot
        xvec = contourvecs(2,:);     % x-axis vector (normalized) for contour plot
        yvec = contourvecs(3,:);     % y-axis vector (normalized) for contour plot

        % Find angle between xvec and yvec by making two new vectors, xvec_orth
        % and yvec_orth that are defined by orthogonal basis vectors
        if or(strcmp(type,'tetragonal') == 1, or(strcmp(type,'orthorhombic') == 1,...
              strcmp(type,'cubic') == 1))
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
        figure(n_mnr+6);

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
    %     disp(xvals_to_delete)
    %     disp(yvals_to_delete)



        % Make contour plot
        cblabel2 = zeros(n_mnr,25); 
        cblabel3 = zeros(n_mnr,6);
        for in = 1:n_mnr
            cblabel2(in,:) = round(linspace(isovalue(in),polmaxa(in),25),3);
            cblabel3(in,:) = round(linspace(isovalue(in),polmaxa(in),6),3);
            figure(n_mnr+6); hold on; 

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
        xlim([0,1])
        ylim([0,1])

    end

    toc
end

% Subfunctions

function [v,c] = get_voronoi(cell_d,angles,phase)
    % get number of atoms 
    atomidx = get_atomloc(phase);
    n_atoms = size(atomidx,1);
    % generate large number of neighbours
    neighbours = gen_neighbours(cell_d,angles,phase);
    % find voronoi cells for all points
    [v,c] = voronoin(neighbours); % v is 
    % reduce to only the cells relevant to the unit cell points, not their
    % neighbours
    %c = c(1:n_atoms);
end

function neighbours = gen_neighbours(cell_d,angles,phase)
    basis = get_basis(cell_d,angles);
    atomidx = get_atomloc(phase);
    neighbour_idx = [1,0,0;0,1,0;0,0,1;1,1,0;1,-1,0;1,0,1;1,0,-1;0,1,1;0,1,-1;1,1,1;-1,1,1;1,-1,1;1,1,-1];
    
    neighboursTemp = atomidx;
    for i = 1:size(neighbour_idx,1)
        neighboursTemp = cat(1,neighboursTemp, atomidx + neighbour_idx(i,:));
        neighboursTemp = cat(1,neighboursTemp, atomidx - neighbour_idx(i,:));    
    end
    neighbours = neighboursTemp; % if you just want to return the idx
    %neighbours = round(((basis')*(neighboursTemp'))',7);
    
end

function basis = get_basis(cell_d,angles)
    
    % construct basis vectors in cartesian coordinates based on lattice
    % parameters and angles. Generalized for triclinic.
    
    % by convention, align the first axis, with length a, with the x axis 
    
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
    
end

function draw_lattice(cell_d,angle,thick,box_clr)
    % Drawing the Unit Cell Outline

    axes = gca;
    x_pos(1) = 0;
    x_pos(2) = cos(angle(3))*cell_d(2);
    x_pos(3) = cell_d(1)+ cos(angle(3))*cell_d(2);
    x_pos(4) = cell_d(1);
    x_pos(5) = 0 + cos(angle(1))*cell_d(3);
    x_pos(6) = cos(angle(3))*cell_d(2) + cos(angle(1))*cell_d(3);
    x_pos(7) = cos(angle(3))*cell_d(2) + cos(angle(1))*cell_d(3) + cell_d(1);
    x_pos(8) = cell_d(1) + cos(angle(1))*cell_d(3);

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
    y_pos(2) = sin(angle(3))*cell_d(2);
    y_pos(3) = sin(angle(3))*cell_d(2);
    y_pos(4) = 0;
    y_pos(5) = 0 + cos(angle(2))*cell_d(3);
    y_pos(6) = sin(angle(3))*cell_d(2) + cos(angle(2))*cell_d(3);
    y_pos(7) = sin(angle(3))*cell_d(2) + cos(angle(2))*cell_d(3);
    y_pos(8) = 0 + cos(angle(2))*cell_d(3);

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
    z_pos(5) = cell_d(3)*sin(angle(1))*sin(angle(2));
    z_pos(6) = cell_d(3)*sin(angle(1))*sin(angle(2));
    z_pos(7) = cell_d(3)*sin(angle(1))*sin(angle(2));
    z_pos(8) = cell_d(3)*sin(angle(1))*sin(angle(2));

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

    line(X1,Y1,Z1,'color',box_clr,'LineStyle','-','LineWidth',thick)

end

function fixed_val = fit_in_cell(val,gridmax)
% exists to correct for any grid coordinate that is outside the bounds of a unit cell.
% Accepts any grid coordinate, and returns the grid coordinate *inside* the unit cell 
% that matches that grid coordinate. For example, if the grid is discretized into 64 steps in a 
% given Cartesian direction (gridmax = 64), and we have a coordinate at -30 
% (val = -30), this function will output a value of 34. If the coordinate is 288, the function will output 32.
    fixed_val = val;
    if fixed_val < 1
        while fixed_val < 1
            fixed_val = fixed_val + gridmax;
        end
    elseif fixed_val > gridmax + 1
        while fixed_val > gridmax + 1
            fixed_val = fixed_val - gridmax;
        end
    end
end