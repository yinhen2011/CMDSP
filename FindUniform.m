% This function is used to determine whether the binary pattern is uniform or not, 
% and if it is, it is set to 1, otherwise it is set to 0
% It is from the Uniformity Optimization Mechanism (UOM) in the paper Enhanced Texture Classification through a Completed Multi-Domain Shrinkage Pattern by Bin Li*, Xiaochun Xu,* , and Q.M.Jonathan Wu,
% which is submitted to The Visual Computer

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