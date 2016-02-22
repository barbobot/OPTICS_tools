function [outputData, aper, p] = makeStdZernikes(Z,n, dataRotate)
%
b = n/2;
v = linspace(-b,b,n+1); v1 = v(1:end-1); v2 = v(2:end);
[x y] = meshgrid(v1,v2); x = x./(b); y = y./(b); 
p = sqrt( x.^2 + y.^2); 
A = atan2(x,y);
A = imrotate(A,180+dataRotate);
data(1,1:n,1:n)   =   	Z(1) * 	ones(n,n)	;
data(2,1:n,1:n)   =   	Z(2) *  4.^(1/2)* (p) .*(cos (A))	;
data(3,1:n,1:n)   =   	Z(3) * 	4.^(1/2)* (p) .*(sin (A))	;
data(4,1:n,1:n)   =   	Z(4) * 	3.^(1/2)* (2*p.^2 - 1)	;
data(5,1:n,1:n)   =   	Z(5) * 	6^(1/2)* (p.^2) .*(sin (2*A))	;
data(6,1:n,1:n)   =   	Z(6) * 	6^(1/2)* (p.^2) .*(cos (2*A))	;
data(7,1:n,1:n)   =   	Z(7) * 	8^(1/2)* (3*p.^3 - 2*p) .*(sin (A))	;
data(8,1:n,1:n)   =   	Z(8) * 	   8.^(1/2)* (3*p.^3 - 2*p) .*(cos (A))	;
data(9,1:n,1:n)   =   	Z(9) * 	   8.^(1/2)* (p.^3) .*(sin (3*A))	;
data(10,1:n,1:n)   =   	Z(10) * 	   8.^(1/2)* (p.^3) .*(cos (3*A))	;
data(11,1:n,1:n)   =   	Z(11) * 	   5.^(1/2)* (6*p.^4 - 6*p.^2 + 1)	;
data(12,1:n,1:n)   =   	Z(12) * 	  10.^(1/2)* (4*p.^4 - 3*p.^2).*(cos (2*A))	;
data(13,1:n,1:n)   =   	Z(13) * 	  10.^(1/2)* (4*p.^4 - 3*p.^2) .*(sin (2*A))	;
data(14,1:n,1:n)   =   	Z(14) * 	  10.^(1/2)* (p.^4) .*(cos (4*A))	;
data(15,1:n,1:n)   =   	Z(15) * 	  10.^(1/2)* (p.^4) .*(sin (4*A))	;
data(16,1:n,1:n)   =   	Z(16) * 	  12.^(1/2)* (10*p.^5 - 12*p.^3 + 3*p) .*(cos (A))	;
data(17,1:n,1:n)   =   	Z(17) * 	  12.^(1/2)* (10*p.^5 - 12*p.^3 + 3*p) .*(sin (A))	;
data(18,1:n,1:n)   =   	Z(18) * 	  12.^(1/2)* (5*p.^5 - 4*p.^3) .*(cos (3*A))	;
data(19,1:n,1:n)   =   	Z(19) * 	  12.^(1/2)* (5*p.^5 - 4*p.^3) .*(sin (3*A))	;
data(20,1:n,1:n)   =   	Z(20) * 	  12.^(1/2)* (p.^5) .*(cos (5*A))	;
data(21,1:n,1:n)   =   	Z(21) * 	  12.^(1/2)* (p.^5) .*(sin (5*A))	;
data(22,1:n,1:n)   =   	Z(22) * 	   7.^(1/2)* (20*p.^6 - 30*p.^4 + 12*p.^2 - 1)	;
data(23,1:n,1:n)   =   	Z(23) * 	  14.^(1/2)* (15*p.^6 - 20*p.^4 + 6*p.^2) .*(sin (2*A))	;
data(24,1:n,1:n)   =   	Z(24) * 	  14.^(1/2)* (15*p.^6 - 20*p.^4 + 6*p.^2) .*(cos (2*A))	;
data(25,1:n,1:n)   =   	Z(25) * 	  14.^(1/2)* (6*p.^6 - 5*p.^4) .*(sin (4*A))	;
data(26,1:n,1:n)   =   	Z(26) * 	  14.^(1/2)* (6*p.^6 - 5*p.^4) .*(cos (4*A))	;
data(27,1:n,1:n)   =   	Z(27) * 	  14.^(1/2)* (p.^6) .*(sin (6*A))	;
data(28,1:n,1:n)   =   	Z(28) * 	  14.^(1/2)* (p.^6) .*(cos (6*A))	;
data(29,1:n,1:n)   =   	Z(29) * 	  16.^(1/2)* (35*p.^7 - 60*p.^5 + 30*p.^3 - 4*p) .*(sin (A))	;
data(30,1:n,1:n)   =   	Z(30) * 	  16.^(1/2)* (35*p.^7 - 60*p.^5 + 30*p.^3 - 4*p) .*(cos (A))	;
data(31,1:n,1:n)   =   	Z(31) * 	  16.^(1/2)* (21*p.^7 - 30*p.^5 + 10*p.^3) .*(sin (3*A))	;
data(32,1:n,1:n)   =   	Z(32) * 	  16.^(1/2)* (21*p.^7 - 30*p.^5 + 10*p.^3) .*(cos (3*A))	;
data(33,1:n,1:n)   =   	Z(33) * 	  16.^(1/2)* (7*p.^7 - 6*p.^5) .*(sin (5*A))	;
data(34,1:n,1:n)   =   	Z(34) * 	  16.^(1/2)* (7*p.^7 - 6*p.^5) .*(cos (5*A))	;
data(35,1:n,1:n)   =   	Z(35) * 	  16.^(1/2)* (p.^7) .*(sin (7*A))	;
data(36,1:n,1:n)   =   	Z(36) * 	  16.^(1/2)* (p.^7) .*(cos (7*A))	;
data(37,1:n,1:n)   =   	Z(37) * 	   9.^(1/2)* (70*p.^8 - 140*p.^6 + 90*p.^4 - 20*p.^2 + 1)	;

aper = p <= 1;
outputData = squeeze(sum(data));
%imagesc(outputData);