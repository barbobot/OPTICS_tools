function [squared] = decodeSquareData(n,data)

v = linspace(-1,1,n);
[x y] = meshgrid(v,v);
x = reshape(x,1,n^2); y = reshape(y,1,n^2);


r = sqrt(x.^2 + y.^2);
rNaN = (r>1);

count = 1
for i=1:length(rNaN);
    if (rNaN(i)==0)
        squared(i) = data(count);
        count = count+1;
    else
        squared(i) = NaN;
    end
end

squared = reshape(squared,n,n);
