function [z_] = afocalConjugateCalculator(z1,f,t)

z = z1;
zn_(1,:) = z*f(1)./(z+f(1));
for i = 2:length(f)
    z(i,:) = -t(i-1)+zn_(i-1,:);
    zn_(i,:) = z(i,:)*f(i)./(z(i,:)+f(i));
end

z_ = zn_(i,:);
figure
plot(z1,z_)

