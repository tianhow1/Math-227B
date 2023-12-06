function [Pz]= Lagrange_interpolation(x,y,z)

n = length(x);

Pz =0;
for i = 1:n
    Pz = Pz +y(i)*Lagrange_polynomial_basis(x,z,i);
end





