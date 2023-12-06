function [Sz]=natural_spline_polynomial(a,b,c,d,x,z)
n=numel(x);
Sz=a(1);
for i=1:n-1
    if((z>x(i))&(z<=x(i+1)))
        Sz=a(i)+b(i)*(z-x(i))+c(i)*(z-x(i))^2+d(i)*(z-x(i))^3;
    end
end
end