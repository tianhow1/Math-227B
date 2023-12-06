function [b] = Divided_difference(x,y)
n = numel(x);
F(:,1) = y; 

for i = 2:n
    for j = 1:i-1
        F(i , j+1) = (F(i,j) - F(i-1,j))/(x(i)-x(i-j));
    end
end

b = diag(F);

end

