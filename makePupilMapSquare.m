function [P_x P_y zeroIndex] = zPupilMapping(n)
v = linspace(-1,1,n);
[x y] = meshgrid(v,v);
x = reshape(x,1,n^2); y = reshape(y,1,n^2);


r = sqrt(x.^2 + y.^2);
P_x = x.*(0./((abs(r)<=1))+1);P_x(isnan(P_x)) = [];
P_y = y.*(0./((abs(r)<=1))+1);P_y(isnan(P_y)) = [];
zeroIndex = (length(P_x)+1)/2;