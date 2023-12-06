function  [Lx]=Lagrange(x,y,j)
[alpha] = Divided_difference(x,y);
Lx = zeros(size(j));
[Lx] = Newton_polynomial(alpha,x,j);
end