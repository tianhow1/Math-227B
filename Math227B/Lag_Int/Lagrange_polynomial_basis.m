function [L_i]= Lagrange_polynomial_basis(x,z,i)

n = length(x);

L_i =1;

for j = 1:n
    if j~=i
        L_i = L_i.*(z-x(j))/(x(i)-x(j));
    end
end



