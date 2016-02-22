function [d d_ P_total] = gaussianReduction(R,n,t, n_medium)
%%
%% Uses Gaussian Reduction to find distance to front and rear
%% principle planes - measured from first and last surface, respectively
%% units are in whatever is input
%% 
%INPUTS - R : radii of curvature (vector of i+1 surfaces)
%%        n : indices or refraction (vector, length i)
%         t : thicknesses between surfaces (vector, length i)
%%        n_medium :  (usually =1, air) index of medium surfaces are surrounded by
%OUTPUTS
%         d : distance from first surface to Principal Plane (object space)
%         d_ : distance from last surface to Rear Principal Plane (image
%         space)
%         P_total : Total power of system (1/focal length)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[p q] = size(n);
if (q==1); n = n';end

%works front to back to reduce surfaces
n = [n_medium, n, n_medium];
P_1 = (n(2)-n(1))/R(1)  %calc power of first surface
t_ = t;  %t_ gets adjusted for principal plane shifts
for i = 1:length(t)
    P_2 = (n(i+2)-n(i+1))/R(i+1)   %calc power of 2nd surface
    P_n = P_1 + P_2 - P_1*P_2*t_(i)/n(i+1);   %Gaussian reduction equation
    d(i) = n(1)*(P_2/P_n)*t_(i)/n(i+1);
    d_(i) = -n(i+2)*(P_1/P_n)*t_(i)/n(i+1);
    if i< length(t)  
        t_(i+1) = t_(i+1) - d_(i);  %adjusts thickness for principal plane shift
    end
    P_1 = P_n
end

effl = 1/P_n
d = sum(d);
d_= d_(end);
P_total = P_n;
    
