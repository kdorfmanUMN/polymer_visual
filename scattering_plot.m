% This script calculates the predicted scattering intensities for the
% structure in the PSCF output by taking the Fourier transform of the
% real-space compositions.

function scattering_plot(R,x,y,z,options)
    
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
        
        % If R is a data array, then the real-space coordinates
        % corresponding to those data must be provided as well. For a 2D
        % system, x and y are needed, while z is also needed in 3D. For a
        % system discretized into an i x j x k grid, these arrays must
        % have size (i,j,k), where x(i,j,k) corresponds to the x-coordinate
        % of the data points in R(i,j,k,:), and so on for y and z.
        
        % x and y are only made optional to allow for the user to run
        % individual_profiles(filename) without causing an error.
        % If the first input parameter is a filename (a string), x, y, and 
        % z are deleted and regenerated by the function read_rgrid.
        % However, if R is a data array and x and y are not provided, an
        % error will occur.
        x = []
        y = []
        z = []

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % The rest of the inputs are optional name-value pair inputs:
        
        % The array scatterers is where the user can specify which
        % species to use as the scattering objects for our Fourier
        % transform. If an array is provided (e.g. [1,2]), the
        % compositions of species 1 and 2 will be added together, and the
        % sum will be used in the Fourier Transform.
        options.scatterers = 1;

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
        
        % writefile is a filename to which the scattering data can be
        % written. This is very useful if you want to do further analyses
        % of the data.
        options.writefile = "";
        
        % hkls is an array in which each row corresponds to one scattering
        % peak to include on the plot (if it is an actual peak). If hkl is
        % not provided, the code below generates this array to include all
        % values of h, k, and l that are 5 or lower. 
        options.hkls = [];

        % fontsize specifies the FontSize parameter for the axis on which
        % data are plotted. Default value is 14.
        options.fontsize = 14;

        % title specifies a string to use as the figure title. Default is
        % no title ("")
        options.title = "";

        % fieldId is an optional index to specify which field to read from 
        % an FTS simulation output file. Default = 0. If the input file is 
        % not an FTS simulation output file, or if R, x, y, and z are 
        % provided as input data arrays, this parameter does nothing.
        options.fieldId = 0;
        
        % theta_plot is a boolean, where we plot our scattering peaks as a
        % function of 2θ if it is true. Otherwise, we plot it as a function
        % of q.
        options.theta_plot = false;
        
        % units is a string, which represents the units of length used in
        % defining x, y, and z. This will be used in the x-axis label of
        % the plot, if provided. So, if units = 'Å', the x-axis label will
        % read 'q [Å^{-1}]'.
        options.units = "";

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
        % plotted are not under a thin film constraint, or have already 
        % been corrected as desired by the user.
        options.film_params;
        
        % Option to turn off the labeling of the peaks
        options.no_labels = false;
        
    end
    
    % Ensure that the code below can access our utilities
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(filepath+"/utilities")
    
    % if a filename is passed to the function, read data from that file
    if ischar(R) || isstring(R) 
        
        clear x y z; % We will determine x, y, and z from the rgrid file
                
        % Read data from file
        [R,x,y,z] = read_rgrid(R,options.fieldId);
        
        % Get lattice basis vectors
        basis = [x(end,1,1),y(end,1,1),z(end,1,1);
                 x(1,end,1),y(1,end,1),z(1,end,1);
                 x(1,1,end),y(1,1,end),z(1,1,end)];
        
    else % get basis and kbasis from x, y, z arrays
        
        basis = [x(end,1,1) y(end,1,1) z(end,1,1);
                 x(1,end,1) y(1,end,1) z(1,end,1);
                 x(1,1,end) y(1,1,end) z(1,1,end)];
    end

    % Apply thin film correction if desired
    if isfield(options,'film_params') && ~isempty(options.film_params)
        [R,x,y,z,basis] = thin_film_correction(R,x,y,z,...
                      options.film_params(1),options.film_params(2),...
                      options.film_params(3),options.film_params(4));
    end

    % Calculate kbasis
    V = abs(dot(basis(1,:),cross(basis(2,:),basis(3,:))));
    kbasis = zeros(3); 
    kbasis(1,:) = round(cross(basis(2,:),basis(3,:))/V,10);
    kbasis(2,:) = round(cross(basis(3,:),basis(1,:))/V,10);
    kbasis(3,:) = round(cross(basis(1,:),basis(2,:))/V,10);


    % Create list of hkl indices to check for scattering peaks if not
    % provided (default list considers all combinations of 0 through 5)
    hkls = options.hkls;
    if isempty(hkls)
        hkls = zeros(6^3,3);
        count = 0;
        for h = -6:6
            for k = -6:6
                for l = -6:6
                    if h ~= 0 || k ~= 0 || l ~= 0
                        count = count + 1;
                        hkls(count,:) = [h,k,l];
                    end
                end
            end
        end
    end

    % Create fgrid, which is a list of coordinates in reciprocal space
    % corresponding to the hkl indices specified in hkls, such that
    % fgrid(i,:) is the vector in reciprocal space described by hkls(i,:)
    % for all i. Also, create q, where q(i) is the magnitude of vector
    % fgrid(i,:).
    fgrid = zeros(size(hkls,1),3);
    q = zeros(size(hkls,1),1);
    for row = 1:size(hkls,1)
        fgrid(row,:) = ((hkls(row,1) * kbasis(1,:)) + ...
                        (hkls(row,2) * kbasis(2,:)) + ...
                        (hkls(row,3) * kbasis(3,:)));
        q(row) = (2*pi) * norm(fgrid(row,:));
    end
    
    % Get a 3D data set D where D(i,j,k) is the sum of the compositions of
    % all scatterers at gridpoint i,j,k. This allows users to specify
    % multiple species as scatterers (i.e. both A and C scatter equally).
    D = zeros(size(x));
    for s = options.scatterers
        D = D + R(:,:,:,s);
    end

    % Reshape x, y, z, and R arrays from 3d to 1d, in order to input the
    % real-space data correctly into MATLAB's built-in FFT function.
    xf = reshape(x(1:end-1,1:end-1,1:end-1),[],1);
    yf = reshape(y(1:end-1,1:end-1,1:end-1),[],1);
    zf = reshape(z(1:end-1,1:end-1,1:end-1),[],1);
    rf = reshape(D(1:end-1,1:end-1,1:end-1),[],1);

    % Perform the FFT using nufftn, a function that performs an FFT on
    % n-dimensional data that is non-uniformly distributed in space. This
    % function is used so that the FFT is correct for all types of unit
    % cells; a uniform FFT only works for cubic data.

    % The output array, Y, is a 1d list of complex numbers, where Y(i) is
    % the Fourier transform of our data evaluated at the coordinates
    % fgrid(i,:) in reciprocal space.
    Y = nufftn(rf,[xf,yf,zf],fgrid);
    Y = Y / numel(xf);
    I = Y .* conj(Y); % Intensities of each scattering peak from FFT data

    % Resort q, I, and hkls in order by q value
    [q,sortInd] = sort(q);
    I = I(sortInd);
    hkls = hkls(sortInd,:);

    % Delete all rows of data with negligible scattering intensities.
    % Negligible is defined as being >10 orders of magnitude less than the
    % highest intensity scattering peak.
    counter = 1;
    I_min = max(I) * 1e-10;
    while counter <= size(I,1)
        if I(counter) < I_min
            I(counter) = [];
            q(counter) = [];
            hkls(counter,:) = [];
        else
            counter = counter + 1;
        end
    end
    
    % Write scattering data to a text file, if writefile is provided
    if options.writefile ~= ""
        writefile = fopen(options.writefile,'w');
        filestring = "%i\t%i\t%i\t%.5f\t%.6e\n";
        fprintf(writefile,"h\tk\tl\tq\tI(q)\n");
        for id = 1:size(I,1)
            fprintf(writefile,filestring,hkls(id,:),q(id),I(id));
        end
        fclose(writefile);
    end

    % Next, additively combine peaks that occur at the same q. Also, create
    % a dictionary called hkls_q, where hkls_q{i} is a matrix where each
    % row represents a family of planes that has a scattering peak at q(i).
    % If two sets of planes {hkl} and {h'k'l'} give scattering peaks at
    % q=q(i), then hkls_q{i} will be the following matrix: [h k l; h' k'
    % l']. A plane (xyz) is a member of the family {hkl} if [x,y,z] is some
    % permutation of [h,k,l]. We need hkls_q in order to ensure that all
    % families of planes with non-negligible scattering are indexed
    % correctly on the plot.
    row = 0;
    hkls_q = cell(size(q,1));
    while row < size(q,1)
        row = row + 1;
        row2 = 1;
        hkls_q{row} = hkls(row,:);
        while row2 <= size(q,1)
            if (row ~= row2) && (abs(q(row)-q(row2)) < 1e-6)
                I(row) = I(row) + I(row2);
                I(row2) = [];
                q(row2) = [];

                ismem = false;
                for id = 1:size(hkls_q{row},1)
                    row_hkl = hkls_q{row}(id,:);
                    all_perms = zeros(48,3);
                    all_perms(1:6,:) = perms(row_hkl);
                    all_perms(7:12,:) = perms(row_hkl.*[-1 1 1]);
                    all_perms(13:18,:) = perms(row_hkl.*[1 -1 1]);
                    all_perms(19:24,:) = perms(row_hkl.*[1 1 -1]);
                    all_perms(25:30,:) = perms(row_hkl.*[-1 -1 1]);
                    all_perms(31:36,:) = perms(row_hkl.*[-1 1 -1]);
                    all_perms(37:42,:) = perms(row_hkl.*[1 -1 -1]);
                    all_perms(43:48,:) = perms(row_hkl.*[-1 -1 -1]);
                    if ismember(hkls(row2,:),all_perms,'rows')
                        ismem = true;
                        hkls_q{row}(id,:) = hkls(row2,:);
                        break
                    end
                end

                if ismem == false
                    hkls_q{row}(end+1,1:3) = hkls(row2,:);
                end

                hkls(row2,:) = [];
            else
                row2 = row2 + 1;
            end
        end
    end

    % Apply permutation to hkls_q to match the thin film correction, if the
    % correction included a rotation (this only affects the labels on the
    % plot, not any part of the numerical calculation)
    if isfield(options,'film_params') && ...
       ~isempty(options.film_params) && options.film_params(4)
        if options.film_params(1) == 0 % normalVec == 0
            perm = [2,3,1];
        elseif options.film_params(1) == 1 % normalVec == 1
            perm = [3,1,2];
        else
            perm = [1,2,3]; % normalVec == 2
        end
        for row = 1:size(hkls,1)
            hkls_q{row} = hkls_q{row}(:,perm);
        end
    end

    % Plot results as I(q) vs q, or I(2θ) vs. θ if theta_plot == true
    figure(); 
    hold on; 
    box on; 
    set(gca,'fontsize',options.fontsize)
    if options.title ~= ""
        title(options.title)
    end

    if options.theta_plot % If we want the x-axis to be 2θ
        twotheta = asin(q / max(q)) * 2 * 180 / pi; 

        for row = 1:size(hkls,1) % Loop over peaks, plot them individually
            plot([twotheta(row),twotheta(row)],[I_min,I(row)],'k',...
                 'linewidth',1)

            % Create label for hkl indices of this peak and put it on plot
            hkl_label = "";
            for id = 1:size(hkls_q{row},1)
                if id == 1
                    hkl_label = strcat(hkl_label,"(",...
                                strjoin(string(hkls_q{row}(id,:)),""),")");
                else
                    hkl_label = strcat(hkl_label,", (",...
                                strjoin(string(hkls_q{row}(id,:)),""),")");
                end
            end
            if options.no_labels == false
                hkl = text(twotheta(row),I(row)*1.05,hkl_label,...
                       'fontsize',options.fontsize*0.9);
                set(hkl,"rotation",90)
            end
        end
        
        xlabel('2\theta [°]')

    else % x-axis is q

        for row = 1:size(hkls,1) % Loop over each peak to plot
            plot([q(row),q(row)],[I_min,I(row)],'k','linewidth',1)

            % Create hkl label for this peak and put it on plot
            hkl_label = "";
            for id = 1:size(hkls_q{row},1)
                if id == 1
                    hkl_label = strcat(hkl_label,"(",...
                                strjoin(string(hkls_q{row}(id,:)),""),")");
                else
                    hkl_label = strcat(hkl_label,", (",...
                                strjoin(string(hkls_q{row}(id,:)),""),")");
                end
            end
            if options.no_labels == false
                hkl = text(q(row),I(row)*1.05,hkl_label,...
                       'fontsize',options.fontsize*0.9);
                set(hkl,"rotation",90)
            end

        end
        
        % Make x-axis label, with or without units
        if options.units ~= ""
            xlabel(strcat('q [',options.units,'^{-1}]'))
        else
            xlabel('q')
        end
        
    end

    % Make plot look good
    ylabel("I(q)")
    set(gca,'yscale','log'); 
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1) pos(2) pos(3)*1.5 pos(4)*1.5]); 
    ylim([min(I)*0.5,max(I)*5]); 
    xlims = xlim; xlim([0,xlims(2)]);
    
    % Save figure if a filename is provided
    if options.savefile ~= ""
        [~,~,ext] = fileparts(options.savefile);
        if (ext == ".fig") || (ext == ".m")
            saveas(gcf,options.savefile);
        else
            exportgraphics(gcf,options.savefile,"resolution",...
                           options.resolution);
        end
    end
    
    hold off;

end
