function gausImagePadded = gausSecDer(im, orders , scale)
%gausSecDer calculates a second gaussian derivative of a given image and scale
%       im       = greyscale image
%       orders   = '[x/y][x/y]'  example: 'xx' for second x derivative
%       scale    = scale(std) of gaussian filter
%
% Stefan Marien, 11/02/2016

    % Build kernel grid
    [X, Y]   = meshgrid( -round(3*scale) : round(3*scale) );

    % Build the gaussian 2nd derivatives filters 
    % (X-Y flip due to matlab coordinates)
    switch orders
        case 'xx'
            gausFilter = 1/(2*pi*scale^4) * (Y.^2/scale^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*scale^2));
        case 'xy'
            gausFilter = 1/(2*pi*scale^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*scale^2));
        case 'yy'
            gausFilter = 1/(2*pi*scale^4) * (X.^2/scale^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*scale^2));

    end

    % Convolution
    gausImageInit = conv2(im, gausFilter, 'valid');
    gausImagePadded = padarray(gausImageInit,(size(im)-size(gausImageInit))./2);
    
end

