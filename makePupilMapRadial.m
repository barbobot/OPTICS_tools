function [P_x P_y zeroIndex] = makePupilMapRadial(arms,rings)

theta = (2*pi/arms)*(0:arms-1) + pi/2;
r = (1/(rings-1)) *(1:rings-1);
%
for i = 1:length(r)
        px(i,:) = r(i)*cos(theta);
        py(i,:) = r(i)*sin(theta);
end
[p q] = size(px);
P_x = reshape(px,1,p*q); P_x = [0 .0001 0    P_x];
P_y = reshape(py,1,p*q); P_y = [0  0   .0001 P_y];
zeroIndex = 1;
%