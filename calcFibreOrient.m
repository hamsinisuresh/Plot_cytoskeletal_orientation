function [ degrees, varargout] = calcFibreOrient( inImage, scales , threshold)
%calcFibreOrient: function to calculate fibre orientation
%
%   Example: degrees = calcFibreOrient( inImage, scales , threshold)
%       return the found orienations in degrees(0-180)
%
%   Example:[ degrees, multiScaleHessian, vessMat, finalVess, dirs, pos] ...
%               = calcFibreOrient( im_green, scales , threshold);
%       return the degrees and all the intermediary results
%
% Stefan Marien, 11/02/2016

    imDims = size(inImage);
    
    %% -- Multiscale Hessian

    % Pre-allocation
    multiScaleHessian = zeros(size(inImage,1), size(inImage,2), size(scales,2), 3 );

    % Loop through scales
    for i = 1:size(scales,2)
        multiScaleHessian(:,:,i,1) = gausSecDer(inImage, 'xx', scales(i) );
        multiScaleHessian(:,:,i,2) = gausSecDer(inImage, 'xy', scales(i) );
        multiScaleHessian(:,:,i,3) = gausSecDer(inImage, 'yy', scales(i) );
    end

    %% -- Pixelwise vesselness calculation by eigenvalue decomposition

    % Pre-allocation
    vessMat = zeros(size(multiScaleHessian, 1), size(multiScaleHessian, 2), size(multiScaleHessian, 3) );
    eigVecs = zeros(size(multiScaleHessian, 1), size(multiScaleHessian, 2), size(multiScaleHessian, 3), 2);

    % Loop through each scale
    for iScale = 1:size(scales,2)

        % Calculate c value for each image (see vesselness formula)
        hessNorms = sqrt( multiScaleHessian(:,:,iScale,1).^2 + multiScaleHessian(:,:,iScale,2).^2 + multiScaleHessian(:,:,iScale,3).^2 ) ;
        cValue = 0.5 * max( abs( hessNorms(:) ) );

        % Loop through image     
        for y = 1:imDims(2)
            for x = 1:imDims(1)

                % Reconstruct Hessian mat (form row to mat)
                curHessMat = [ multiScaleHessian(x,y,iScale,1), multiScaleHessian(x,y,iScale,2);
                               multiScaleHessian(x,y,iScale,2), multiScaleHessian(x,y,iScale,3) ];

                % Calculate frangi vesselness value
                [vess, eigVec] = calcFrangiVesselness(curHessMat, cValue);
                vessMat(x,y,iScale) = vess;

                % Save largest eigenvectors (principal direction)
                eigVecs(x,y,iScale,:) = eigVec;

            end
        end
    end
    
    %% -- Select maximum vesselness across scales  
    bestScales = zeros(size(multiScaleHessian,1),size(multiScaleHessian,2));

    for y = 1:imDims(2)
        for x = 1:imDims(1)

            % Find scale with highest vesselness, 'last' makes it take the
            % highest scale when several scales have the same vesselness
            maxVes = max(vessMat(x,y,:));
            bestScales(x,y) = find(vessMat(x,y,:) == maxVes, 1, 'last');

        end
    end

    % Build vesselness image
    finalVess = max( vessMat, [], 3) ;

    % Threshold vesselness image 
    finalVessThres = im2bw(finalVess, threshold) .* finalVess;
    
    %% -- Extract directions after threshold

    % Pre-allocation
    dirs = zeros(2, sum(finalVessThres(:) > 0.995) );
    pos  = zeros(2, sum(finalVessThres(:) > 0.995) );
    counter = 1;

    % Loop through image
    for y = 1:imDims(2)
        for x = 1:imDims(1)

            % Use only pixels from thresholded vesselness image
            if finalVessThres(x,y) > 0

                % Get best scale
                iScale = bestScales(x,y);

                % Export to list variables
                dirs(:,counter) = [eigVecs(x,y,iScale,2) ; eigVecs(x,y,iScale,1)];
                pos(:,counter) = [x ; y];

                counter = counter + 1;
            end

        end
    end

    %% -- Orientation statistics

    % Pre-allocation
    degrees =  zeros(1, size(pos,2));

    % Loop though each vessel pixel
    for i = 1:size(pos,2)

        % Calculate angle
        degrees(i) = atan2d( dirs(1,i), dirs(2,i) );

        % Combine into 180 degree range
        if degrees(i) >= 180
            degrees(i) = degrees(i) - 180;

        elseif degrees(i) <= 0
            degrees(i) = degrees(i) + 180;        

        end

    end
    
    % Handle optional out arguments
    if nargout~=1
        varargout{1} = multiScaleHessian ;
        varargout{2} = vessMat ;
        varargout{3} = finalVess ;
        varargout{4} = dirs ;
        varargout{5} = pos ;
        varargout{6} = degrees ;
    end
    
end

