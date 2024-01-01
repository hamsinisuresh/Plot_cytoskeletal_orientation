function [vesselness,eigVec] = calcFrangiVesselness( hessMat, cValue)
%calcFrangiVesselnss, calculates the frangi vesselness from a hessian
%matrix
%
% Stefan Marien, 11/02/2016
    
        % Calculate eigensystem
        [eigVecs, eigVals]= eig(round(hessMat,4,'significant'));
        eigVec = eigVecs(:,1);
        eigVals = [eigVals(1,1), eigVals(2,2)];
        
        l1 = eigVals(2);
        l2 = eigVals(1);
        
        % Calculate vesselness
        if l2 >=0
            vesselness = 0;
        else
            vesselness = exp(-2*( l1 / l2 )^2 ) * ...
                (1 - exp( -(l1^2 + l2^2) / 2*((cValue^2)) ) );
        end

end

