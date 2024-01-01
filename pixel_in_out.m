%% function to check if a point lies inside triangle
%
% Hamsini Suresh, 27/04/2018
%  
%{
    % for faster execution of this function, produce the mex version by running the following commands at prompt:

    double1Da = coder.typeof(double(1), [inf, 1], [false, false]);
    double2Da = coder.typeof(double(1), [inf, 2], [false, false]);
    double2Db = coder.typeof(double(1), [inf, 3], [false, false]);    
    sArgs = {double(1), double(1), double1Da, double2Da, double2Db};
    codegen pixel_in_out -args sArgs

%}
function [pixel_Tria] = pixel_in_out(pixel_rows, pixel_cols, angle2_col, Nodes, Tria)
    nTria = size(Tria,1);
    pixel_Tria = nan(pixel_rows, pixel_cols);
    angle2 = reshape(angle2_col, pixel_rows, pixel_cols);
    
    for ii = 1:pixel_cols
        for jj = 1:pixel_rows  
            if ~isnan(angle2(jj,ii))  % pixel exists
                Cnuc = [ii, -jj];
                % put each pixel into a mesh element or nan
                for kk = 1:nTria
                    tri1 = Nodes(Tria(kk,1),:);
                    tri2 = Nodes(Tria(kk,2),:);
                    tri3 = Nodes(Tria(kk,3),:);    

                    % mid-points
                    tri12 = tri2 - tri1;
                    tri23 = tri3 - tri2;
                    tri31 = tri1 - tri3;

                    n12 = [tri12(2), -tri12(1)];
                    if dot(n12, -tri31) > 0
                        n12 = -n12;
                    end
                    n23 = [tri23(2), -tri23(1)];
                    if dot(n23, -tri12) > 0
                        n23 = -n23;
                    end
                    n31 = [tri31(2), -tri31(1)];
                    if dot(n31, -tri23) > 0
                        n31 = -n31;
                    end

                    if (dot(n12, (Cnuc-tri1)) < 0) && (dot(n23, (Cnuc-tri2)) < 0) && (dot(n31, (Cnuc-tri3)) < 0) % pixel inside this element
                        pixel_Tria(jj,ii) = kk;    
                    end
                    
                end
            end
        end
    end   
end