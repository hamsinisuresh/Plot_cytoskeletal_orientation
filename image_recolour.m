function image_recolour(start_id)
%% colour pixels within mesh by local 'average' orientation
%
% Hamsini Suresh, 27/04/2018
%
clc

strpaths = {'SLB1','SLB2','SLB3','SLC1','SLC2',...
    'SLC3','SLC4','SLE1','SLE2','Homog'};
numcells = [45; 60; 48; 48; 57; 48; 50; 61; 57; 52];

for ru = start_id:start_id
str2 = strcat('Plot_cytoskeletal_orientation\raw_images\',strpaths(ru));
str3 = char(str2);
addpath(str3);
cellnum1 = numcells(ru);

opcsk2 = zeros(cellnum1,1);

for ty = 1:cellnum1
    fprintf('Cell %d\n', ty);    
    
    strfile = sprintf('%s\\csk_%d.mat', str3, ty);  % file with coordinates of mesh covering cell
    load(strfile, 'smooth_xy', 'Nodes', 'Tria');
    
    strfile2 = sprintf('%s\\cskdiff_angle_%d.mat', str3, ty); % file with SF orientations
    load(strfile2, 'angle2', 'angle3', 'binnum');   % angle2 has angles in finer bins
    pixel_rows = size(angle2,1);
    pixel_cols = size(angle2,2);
    nTria = size(Tria,1);
    
    angle2_col = angle2(:);
    angle3_col = angle3(:);
    avgcsk = 112.2605;
    
	% use mex version to run this function faster
    [pixel_Tria] = pixel_in_out(pixel_rows, pixel_cols, angle2_col, Nodes, Tria);
    
    % ----- In which element does the pixel lie? ----- %
    % point inside triangle                            
    disp('Sorted into elements')   
    
    % find all pixels in a given element
    imp_angle2 = zeros(pixel_rows, pixel_cols);
    imp_angle3 = zeros(pixel_rows, pixel_cols);
    phi_Tria = nan(nTria,1);   % -90 to +90
    color_Tria = nan(nTria,3);   % -90 to +90
    
    hsv_size = 4096;  
    myColorMap2 = hsv(hsv_size); 
    myColorMap3 = hsv(binnum);
%     myColorMap(1,:) = 0;

    phi2 = linspace(-90,90,hsv_size)';       
    pixel_Tria_col = pixel_Tria(:);    
    pixel_Tria_col2 = pixel_Tria(:);

    for rr = 1:nTria
        B = find(pixel_Tria_col == rr);   % all pixels in a given element            
        if ~isempty(B)   % if atleast one pixel exists in element            
            imp_angle2 = angle2_col(B);   % all angles2 within one element            
            phi_Tria(rr) = mean(imp_angle2);   % mean on finer bins
            pixel_Tria_col2(B) = mean(imp_angle2);
        end
        
        if isnan(phi_Tria(rr))    % no pixels => white
            color_Tria(rr,:) = myColorMap2(hsv_size/2,:);                      
        else            
            [~,L] = min(abs(phi_Tria(rr)-phi2)); % csk diff
            color_Tria(rr,:) = myColorMap2(L,:);   
        end          
    end    
    pixel_angles_re = reshape(pixel_Tria_col2, pixel_rows, pixel_cols);    

    binnum = 36;    
    phi3 = linspace(-90,90,binnum)';    
    fig1 = ones(pixel_rows, pixel_cols, 3);
    fig2 = ones(pixel_rows, pixel_cols, 3);
    fig3 = ones(pixel_rows, pixel_cols, 3);
        
    for ii=1:pixel_rows
        for jj=1:pixel_cols 
            % Original figure
            if isnan(angle2(ii,jj))
                fig1(ii,jj,:) = [1 1 1];  % [0 0 0]            
            else   
                [~,L] = min(abs(angle2(ii,jj)-phi2-90+avgcsk)); % csk diff
                fig1(ii,jj,:) = myColorMap2(L,:);                                        
            end

            % meshed coloured figure
            if isnan(pixel_angles_re(ii,jj))
                fig2(ii,jj,:) = [1 1 1];  % [0 0 0]            
            else   
                [~,L] = min(abs(pixel_angles_re(ii,jj)-phi2-90+avgcsk)); % csk diff
                fig2(ii,jj,:) = myColorMap2(L,:);                                        
            end
        end
    end


    % find all pixels which differ by more than 30 deg
    new_angles = abs(pixel_angles_re - angle2);
    final_angles = angle2;

    ctr = 0;
    ctr1 = 0;
    for ii=1:pixel_rows
        for jj=1:pixel_cols 
            % meshed RECOLOURED figure
            if isnan(angle2(ii,jj))
                fig3(ii,jj,:) = [1 1 1];  % [0 0 0]            
            elseif new_angles(ii,jj) >= 50 
                ctr = ctr + 1;  % how many speckles?
                final_angles(ii,jj)= pixel_angles_re(ii,jj);
                fig3(ii,jj,:) = fig2(ii,jj,:);  
            else
                ctr1 = ctr1 + 1;
                fig3(ii,jj,:) = fig1(ii,jj,:);
            end
        end
    end
    
        final_angles1 = final_angles(:);
        final_angles2 = final_angles1(~isnan(final_angles1));

    opcsk2(ty) = sqrt( (mean(cosd(2*final_angles2)))^2 + (mean(sind(2*final_angles2)))^2 );
    fprintf('%d  op-csk * log-magnitude: %.4f\n', ty, opcsk2(ty));

    % SCALE BAR
    rgbImage3 = ones(pixel_rows+50, pixel_cols+200, 3);
    rgbImage3(51:pixel_rows+50, 51:pixel_cols+50, :) = fig3;
%     rgbImage3(20, 5:117, :) = 0;  % scale bar for 20x; 113 px is 60 micron
    rgbImage3(20, 5:230, :) = 0;  % scale bar for 40x; 226 px is 60 micron
    
    
% --- Visualize image with stress-fibres recolored by "adjusted" orientations
    figure(4); clf
    set(gcf,'color','w'); 
    imshow(rgbImage3);

    pause(2)
    fprintf('done\n\n')

    sFiles = sprintf('image_adjusted_csk_%d.svg', ty);  % tif or svg
    saveas(gca, sFiles);
end
end
end