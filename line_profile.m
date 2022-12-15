% This function creates a figure that represents the composition of each
% species in a polymer system as a function of position along a line
% through the unit cell. 

% It is generalized so that the endpoints of the line (startloc and
% startloc+direc) do not need to fall on a gridpoint, you can place them
% wherever you want. If there are no gridpoints on the line, it will just
% give you an empty plot.

function line_profile(R,direc,startloc,options)

    arguments
        
        % The first parameter of this function, R, is "overloaded". This
        % means that the function accepts either a string or a data array 
        % for the first parameter, and it will work correctly either way. 
        
        % If the first parameter is a string, we assume it is a filename
        % that contains the data we wish to plot. Thus, we will use
        % read_rgrid.m to collect the data that will be plotted.
        
        % If the first parameter is instead a data array, it must contain
        % all of the composition data stored in the rgrid file. For an
        % N-dimensional system, R must be an (N+1)-dimensional array. In
        % 3D, R(i,j,k,l) gives the composition of species l at gridpoint
        % (i,j,k). 
        R
        
        % The other required parameter for this function is direc, which is
        % the direction of the vector along which we trace the line
        % profile. This is defined in reduced coordinates. The length of
        % the vector direc also defines the length of the line along which
        % our profile is traced.
        direc {mustBeNumeric}
        
        % startloc is an optional third input that specifies the starting
        % point of our line profile. The default position is [0 0 0].
        startloc = [0 0 0]
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % The rest of the inputs are optional name-value pair inputs:

        % savefile is a filename to which the figures will be saved.
        % The file extension provided (e.g. ".fig" or ".png") will be used
        % to determine the type of file to save. Since this function
        % typically generates more than 1 figure, we append mono_label(i)
        % to the end of the filename when saving the composition profile
        % for species i. If left empty (which is the default behavior), the
        % figures are not saved.
        options.savefile = "";

        % resolution is a number that specifies the resolution of the
        % figure that is saved (if options.savefile is specified), in dots
        % per inch (dpi). Default value is 300. If set to 0, file is saved
        % at screen resolution.
        options.resolution = 300;

        % fontsize specifies the FontSize parameter for the axis on which
        % data are plotted. Default value is 10.
        options.fontsize = 14;
        
        % mono_label contains labels for each monomer species. If not
        % specified, we use ["A","B","C",...] as the default behavior.
        options.mono_label = ["A","B","C","D","E","F","G","H"];
                
        % colors defines the color of each line. colors[i,:] is an RGB
        % triplet that is the color for the line profile of species i.
        % The default values defined below correspond to the default colors
        % given in utilities/get_colormaps.
        options.colors = [0,   0.7, 0.9;   %blue
                          0.9, 0,   0;     %red
                          0,   0.9, 0.2;   %green
                          1,   1,   0;     %yellow
                          0.5, 0,   1;     %purple
                          1,   0,   1;     %pink
                          1,   0.5, 0;     %orange
                          0.75,0.75,0.75]; %grey

        % If your SCFT result is a thin film, you should include
        % film_params as an input to apply a thin film correction.
        %
        % film_params is an array with 4 entries. The first 3 entries
        % correspond to the 3 required parameters in pscfpp that are needed
        % to define a Wall object: normalVec, interfaceThickness, and
        % wallThickness. See pscfpp documentation for details about what
        % each of these three parameters means. The fourth entry is a
        % boolean (i.e. 0 for false, 1 for true) that indicates whether or
        % not to rotate the figure to make the z axis orthogonal to the
        % wall. If this film_params input is included, the code will apply
        % a correction to the plot to make the figure look good as a thin
        % film. If it is not included, it is assumed that the data being
        % plotted are not under a thin film constraint.
        options.film_params;
        
        % You must provide coords as an input if the following are true:
        %   1) R, the first input to this function, is a data array (as
        %      opposed to a filename)
        %   2) You wish to apply a thin film correction
        % Coords is merely a cell array containing the x, y, and z arrays
        % that are produced by read_rgrid.m (so, coords = {x, y, z}). These
        % coordinates are required inputs for thin_film_correction.m, so we
        % must have them available if we wish to apply the correction.
        options.coords;
        
    end
    
    % Ensure that the code below can access our utilities
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(filepath+"/utilities")

    % if a filename is passed to the function, read data from that file
    if ischar(R) || isstring(R) 
        
        close all; % close other figures
                
        % Read data from file
        options.coords = cell(1,3);
        [R,options.coords{1},options.coords{2},...
                                options.coords{3}] = read_rgrid(R);
    end
    
    % define dim, grid, and n_mnr
    dim = length(direc); % Note: this does not necessarily need to match
                         % the dim of the system. R always contains 3D data
    grid = ones(1,3);
    for d = 1:dim
        grid(d) = size(R,d) - 1;
    end
    n_mnr = size(R,4); % Note: R contains 3D data, even if dim < 3

    % Make direc and startloc into 3D vectors if not already
    while length(direc) < 3; direc(end+1) = 0; end %#ok<AGROW> 
    while length(startloc) < 3; startloc(end+1) = 0; end %#ok<AGROW> 

    % Apply thin film correction to data if desired
    startloc_og = startloc; % input startloc, before any corrections
    direc_og = direc; % input direc, before any corrections
    if isfield(options,'film_params') && ~isempty(options.film_params)

        if ~isfield(options,"coords") % Make sure we have coords
            error("coords is a required input for this condition")
        end

        R = thin_film_correction(R,options.coords{1},options.coords{2},...
                          options.coords{3},options.film_params(1),...
                          options.film_params(2),options.film_params(3),...
                          false,true);
        
        if options.film_params(4) % If a rotation is desired, then:
            % Change direc and startloc to match the rotation. R is not
            % altered by the rotation in thin_film_correction (only x,
            % y, and z) so we must do this manually in this function
            if options.film_params(1) == 0 % normalVec == 0
                direc = direc(:,[3,1,2]);
                startloc = startloc(:,[3,1,2]);
            elseif options.film_params(1) == 1 % normalvec == 1
                direc = direc(:,[2,3,1]);
                startloc = startloc(:,[2,3,1]);
            end
        end
    end

    % Define variables that characterize the line
    start_coord = startloc .* grid + 1;
    end_coord = start_coord + (direc .* grid);
    dir_vec = end_coord-start_coord;
    
    [~,ind] = max(abs(dir_vec));
    start_data = start_coord;
    end_data = end_coord;
    
    % if start_coord is on a gridpoint
    if max(abs(start_coord-round(start_coord))) < 1e-5 
        x_min = 0;
    else % Find next gridpoint along the line that is before start_coord
        diff = start_data(ind) - floor(start_data(ind));
        start_data = start_data - (diff * direc / direc(ind));
        
        counter = 0;
        
        while max(abs(start_data-round(start_data))) > 1e-5 && ...
              counter < grid(ind)
          
            start_data = start_data - (direc / direc(ind));
            counter = counter + 1;
            
        end
        
        x_min = 1 - (norm(end_coord - start_data) / norm(dir_vec));
    end
    
    % if end_coord is on a gridpoint
    if max(abs(end_coord-round(end_coord))) < 1e-5 
        x_max = 1;
    else % Find next gridpoint along the line that is past end_coord
        diff = ceil(end_data(ind)) - end_data(ind);
        end_data = end_data + (diff * direc / direc(ind));
        
        counter = 0;
        
        while max(abs(end_data-round(end_data))) > 1e-5 && ...
              counter < grid(ind)
          
            end_data = end_data + (direc / direc(ind));
            counter = counter + 1;
            
        end
        
        x_max = norm(end_data - start_coord) / norm(dir_vec);
    end
    
    data_vec = end_data - start_data;
    n_steps = max(abs(data_vec));
    x_vals = linspace(x_min,x_max,n_steps+1);
    gridpoints = zeros(n_steps+1,dim);
    for d = 1:dim
        gridpoints(:,d) = ((0:n_steps)/n_steps) * data_vec(d)...
                          + start_data(d);
    end        

    % Delete rows in gridpoints array that are not on an actual gridpoint
    pt = 1;
    while pt <= size(gridpoints,1) % Loop over rows
        % If gridpoints(pt,:) not on an actual gridpoint
        if max(abs(gridpoints(pt,:)-round(gridpoints(pt,:)))) > 1e-5
            % Delete row
            gridpoints(pt,:) = [];
            x_vals(pt) = [];
        else
            % Make sure gridpoint is inside unit cell, where we have data
            for d = 1:dim
                gridpoints(pt,d) = fit_in_cell(gridpoints(pt,d),grid(d));
            end
            
            % Proceed to next row
            pt = pt + 1;
        end
    end
    gridpoints = round(gridpoints);

    % Make gridpoints into a 3-column array even if dim < 3, to index R
    % properly
    if dim == 2
        gridpoints(:,3) = 1;
    elseif dim == 1
        gridpoints(:,2:3) = 1;
    end
    
    % Get composition data at all gridpoints along line
    line_plot = zeros(n_mnr,length(x_vals));
    for in = 1:n_mnr
        for pt = 1:size(gridpoints,1)
            line_plot(in,pt) = R(gridpoints(pt,1),gridpoints(pt,2),...
                                 gridpoints(pt,3),in);
        end
    end
    
    % Construct plot
    figure(); 
    hold on
    box on
    set(gca,'fontsize',options.fontsize);
    linemarkers = {'s','^','o','p','+','*','d'};
    legend_labels = strings(1,n_mnr);
    
    for in = 1:n_mnr
        plot (x_vals,line_plot(in,:),'color',...
              options.colors(in,:),'marker',linemarkers{in},...
              'markerfacecolor',options.colors(in,:),'markersize',5)
        hold on
        legend_labels(in) = strcat(options.mono_label(in),' block');

    end
    
    if isequal(startloc,[0 0 0])
        title1 = sprintf('Density Profile Along [%g %g %g]',direc_og);
    else
        title1 = sprintf('Density Profile from [%g %g %g] to [%g %g %g]'...
                         ,startloc_og,startloc_og+direc_og);
    end
    title(title1)
    xlabel('r/r_{max}')
    ylabel('\phi(r)')
    xlim([0 1])

    legend(legend_labels)
    legend('Location','northeastoutside')

    % Save figure if a filename is provided
    if options.savefile ~= ""
        [~,~,ext] = fileparts(options.savefile);
        if (ext == ".fig") || (ext == ".m")
            saveas(gcf,options.savefile);
        else
            if ext == ".jpg"
                format = "-djpeg";
            elseif ext == ".tif"
                format = "-dtiff";
            else
                format = strcat("-d", extractAfter(ext,1));
            end
            res = strcat("-r",num2str(options.resolution));
            print(gcf,options.savefile,format,res);
        end
    end
    
    hold off

end
