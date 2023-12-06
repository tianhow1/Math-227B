function  [Lx]=Lagrange_polynomial_interpolation(a,b,f,n,j)
% a is lower bound, b is upper bound, n is the number of point that want to
% interpolate
x = (0:2/(n-1):b-a)+a; %Divide in terval [a,b] in to n-1 area, which generate n point.
y=f(x);%Get the function value
[alpha] = Divided_difference(x,y);
Lx = zeros(size(j));

[Lx] = Newton_polynomial(alpha,x,j);

end