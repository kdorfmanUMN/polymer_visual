
function [map_store,outcolor] = get_colormaps(n)
    
    arguments
        n = 1000 % number of distinct colors in each colormap
    end

    % Drawing the Color Maps
    ncolor = 8;

    % low = low fraction, i.e. light
    color_low = zeros(ncolor,3);
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

    map_store = cell(ncolor,1);

    for in = 1:ncolor
        temp_map = zeros(n,3);
        temp_map(:,1) = linspace(colorpad(in,1,1),colorpad(in,1,2),n); %Red
        temp_map(:,2) = linspace(colorpad(in,2,1),colorpad(in,2,2),n); %Green
        temp_map(:,3) = linspace(colorpad(in,3,1),colorpad(in,3,2),n); %Blue
        map_store{in}=temp_map;
    end
    
end