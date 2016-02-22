function [z_ m] = thinLensEquation(f,z)

z_inv=1/f+1/z;
z_ = 1/z_inv;
m = z_/z 
