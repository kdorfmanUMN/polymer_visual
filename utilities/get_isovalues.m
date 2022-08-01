% Computing the isovalues

function isovalue = get_isovalues(R,dim,options)

    arguments
        R
        dim
        options.plot = false;
        options.colors = [0,   0.7, 0.9;   %blue
                          0.9, 0,   0;     %red
                          0,   0.9, 0.2;   %green
                          1,   1,   0;     %yellow
                          0.5, 0,   1;     %purple
                          1,   0,   1;     %pink
                          1,   0.5, 0;     %orange
                          0.75,0.75,0.75]; %grey
        options.fontsize = 14;
    end

    grid = size(R,1:3) - 1; % Note: grid is always 3D even for dim < 3, as 
                            % defined by the function read_rgrid.m
    n_mnr = size(R,4);
    
    makeplot = options.plot;
    colors = options.colors;
    fontsize = options.fontsize;
    clear options;
    
    weight = ones(1,n_mnr);
    gap_close = 1; %Toggle the white space closing feature (Default = 1)
    plot_grid = zeros(1,dim);
    for ig = 1:dim
        plot_grid(ig) = grid(ig)+1;
    end

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
    
    % Plotting the rescaled volume fractions and isovalues if desired
    if makeplot
        figure(); 
        hold on;
        set(gca,'fontsize',fontsize)
        
        linestyles(1) = {'-.'};
        linestyles(2) = {':'};
        linestyles(3) = {'--'};
        for in = 1:n_mnr
            plot(x_plot,line_new(in,:),'color',colors(in,:))
        end

        clear yy
        yy = ones(n_mnr,length(x_plot));
        for in = 1:n_mnr
            yy(in,:) = yy(in,:)*isovalue_s(in);
        end

        for in = 1:n_mnr
            plot(x_plot,yy(in,:),'color',colors(in,:),...
                 'LineWidth',1+n_mnr-in,'LineStyle',...
                 cell2mat(linestyles(mod(in,3)+1)))
        end

        xlabel('Line Index')
        ylabel('Rescaled Volume Fraction')

        legend_labels1 = strings(1,n_mnr);
        legend_labels2 = strings(1,n_mnr);
        mono_label = char(1,n_mnr);
        for in=1:n_mnr
            mono_label(in) = char(in+'A'-1);
            legend_labels1(in) = strcat(mono_label(in),' block Rescaled Density');
            legend_labels2(in) = strcat(mono_label(in),' block Rescaled Isovalue');
        end


        legend([legend_labels1 legend_labels2])
        legend('Location','best')
        
        hold off;
    end

end