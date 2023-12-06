function [PX] = Newton_polynomial(b,x,X)
n = numel(x);
S = ones(size(X));
PX = b(1)*S;

for i = 2:n
        S = S.*(X-x(i-1));
        PX = PX+ b(i).*S;
end


end

