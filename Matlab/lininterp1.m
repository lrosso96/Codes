function yi = lininterp1(x,Y,xi)
%lininterp1: linearly interpolate function Y defined on X on points Xi.
%
% x must be a strictly increasing vector;
% Y must be a column vector or a matrix with length(X) rows;
% xi is a vector or a scalar;
% yi will be length(xi)-by-size(Y,2);
%
%--------------------------------------------------------------------------
 
F = griddedInterpolant(x,Y);
yi = F(xi);

end