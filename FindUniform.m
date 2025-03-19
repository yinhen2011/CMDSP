function [pattern_i] = FindUniform(pattern,P)
 [x1,y1] = size(pattern);
    pattern_i = zeros(x1,y1);
   
    for j = 1:x1
        for k = 1:y1
            if pattern(j,k) <= P
                pattern_i(j,k) = 1 ;
            else
                pattern_i(j,k) = 0 ;
            end
        end
    end
end