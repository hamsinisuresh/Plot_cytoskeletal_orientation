%% Fibre Orientation Analysis script
%
% Script to analyse fibre orientations in microscopy images
%
% Stefan Marien, 11/02/2016; adapted by Hamsini Suresh, 01/02/2018
%

clear all
close all

%% -- Input Parameters

% Calculation
scales          = [1, 1, 1];  %[0.8, 1, 1.2]; % 
threshold       = 0.6; 
showEachXArrow  = 1;

% Visualisation (1= plot / 0= not plot)
visAngledHisto      = 1 ;  % Angled histogram of orientations
visNormalHisto      = 1 ;  % Normal histogram of orientations
visGausDerImages    = 1 ;  % Gaussian derivative iamges for each scale
visScaleVesselness  = 1 ;  % Vesselness images for each scale
visFinalVesselness  = 1 ;  % Final vesselness images
visOverlayImage     = 1 ;  % Orientations overlayed on input image
HeatOverlayImage     = 1 ;% Heatmap
HeatOverlayImage1 = 1;  %Heatmap

% Add function folder to path
[curDir, ~, ~] = fileparts( mfilename('fullpath') );
addpath( fullfile(curDir, 'functions') );

%% -- Importing

tic

% Import image
im = imread('image_eLoG.tif');

% Take green channel 
im_green = im(:,:,1);

% Convert to integers (ranging from 0 to 255)
im_green = im2double(im_green) .* 255;


%% -- Calculation

% Main calculation:
[ degrees, multiScaleHessian, vessMat, finalVess, dirs, pos] ...
        = calcFibreOrient( im_green, scales , threshold);
toc

%% -- Visualisations

% -- Show angle histogram --
if visAngledHisto == 1
    
    figure;
    
    % Show angle histogram 
    [T1, R1] = rose(degrees./180*pi,100);   %returns the vectors T and R such that 
    R1 = R1./2./sum(R1);                    % Normalise

    rosePlot1 = polar([T1,T1+pi], [R1,R1]);                % is the histogram

    title('Orientation Histogram')
end

% -- Show normal histogram --
if visNormalHisto == 1
    figure;
    
    histogram(degrees,100,'Normalization','probability')
    xlabel('Degrees')
    ylabel('Probability')
    title('Orientation Histogram')

    viewRange = axis;
    axis([0 180 viewRange(3) viewRange(4) ])
    
end

% -- Show Guassian derivatives result --
if visGausDerImages == 1
    figure;
    hold off

    for i = 1:size(scales,2)
        % XX
        subplot( size(scales,2), 3, 1+(i-1)*3 )
        imshow( .5 + multiScaleHessian(:,:,i,1) ./ (2*max(max(abs(multiScaleHessian(:,:,i,1))))) )
        title( strcat(['Scale: ', num2str(scales(i)), ', I_{xx}']) )

        % YY
        subplot( size(scales,2), 3, 2+(i-1)*3 )
        imshow( .5 + multiScaleHessian(:,:,i,2) ./ (2*max(max(abs(multiScaleHessian(:,:,i,2))))) )
        title( strcat(['Scale: ', num2str(scales(i)), ', I_{xy}']) )

        % XY
        subplot( size(scales,2), 3, 3+(i-1)*3 )
        imshow( .5 + multiScaleHessian(:,:,i,3) ./ (2*max(max(abs(multiScaleHessian(:,:,i,3))))) )
        title( strcat(['Scale: ', num2str(scales(i)), ', I_{yy}']) )

    end
end

% -- Visualise vesselness for all scales --
if visScaleVesselness == 1
    figure;
    title('Vesselness')
    for iScale = 1:size(scales,2)
        subplot(1, 3, iScale )
        imshow(vessMat(:,:,iScale))
        title( strcat(['Vesselness, scale: ', num2str(scales(iScale)) ]) )
    end
end

% -- Visualise final vesselness image --
if visFinalVesselness == 1
	figure; imshow(finalVess)
    title('Final vesselness')
end

% -- Visualise overlay --
if visOverlayImage == 1
    figure; 
    imshow(im)
    hold on
    
    quivPlot= quiver( pos(2,1:showEachXArrow:end),...
            pos(1,1:showEachXArrow:end),...
            -dirs(2,1:showEachXArrow:end),...
             dirs(1,1:showEachXArrow:end), ...
            1 );

    quivPlot.Color = 'y';
    
    hold off
end

% -- Visualise heatmap --
hsv_size = 4096;
if HeatOverlayImage == 1
   
    figure; 

posnew=flipud(pos);
xyz=[posnew; degrees];
XYZ=xyz';
Z=accumarray(XYZ(:,[2 1]), XYZ(:,3));
imagesc(Z);
myColorMap = flipud(hsv(hsv_size)); % 4096
myColorMap(1,:) = 0;
colormap(myColorMap);
colorbar


phi2 = linspace(0,180,hsv_size)';
rows1 = size(Z,1);
cols1 = size(Z,2);
fig2 = ones(rows1,cols1,3);
for ii=1:rows1
    for jj=1:cols1      
        if (Z(ii,jj)==0)
            fig2(ii,jj,:) = [0 0 0];
        else
            [~,L] = min(abs(Z(ii,jj)-phi2));
            fig2(ii,jj,:) = myColorMap(L,:); 
        end
    end
end

end


% -- Visualise heatmap --
 if HeatOverlayImage1 == 1
%     
   figure;  
x = XYZ(:,1);
y = XYZ(:,2)
z = XYZ(:,3)   
n=256;
[X, Y] = meshgrid(min(x): max(x), min(y): max(y));
Zint = griddata(x,y,z,X,Y,'cubic');

%// Remove the NaNs for imshow:
Zint(isnan(Zint)) = 0;
imagesc(min(x): max(x),min(y): max(y),Zint)
colormap(myColorMap)
 end 
 
 opcsk = sqrt(mean(cosd(2*z)).^2+mean(sind(2*z)).^2);
 fprintf('Cytoskeletal order parameter: %.4f\n', opcsk);
 
 allmag1 = zeros(length(z),1);
 pos1 = pos';
 for gg=1:length(z)
     allmag1(gg) = im(pos1(gg,1),pos1(gg,2));
 end
csk_order_parameter = sqrt( (sum(allmag1.*cosd(2*z))/sum(allmag1))^2 + (sum(allmag1.*sind(2*z))/sum(allmag1))^2 );
fprintf('csk orientation * magnitude: %.1f\n', csk_order_parameter);
