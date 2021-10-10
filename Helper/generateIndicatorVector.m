function [e] = generateIndicatorVector(index,dimension)
% This function generates an indicator vector e of size (dimension) with a 1 in
% the (index) location and 0's elsewhere
e = zeros(dimension,1);
e(index) = 1;
end

