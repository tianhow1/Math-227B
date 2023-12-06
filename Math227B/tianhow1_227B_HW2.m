clear

close all
%% Question 1
%(a), the answer is in the Divided_difference.m in Lag_Int folder, also in
%the Lag_Int package.
%(b), the answer is in the Newton_polynomial.m in Lag_Int folder, also in
%the Lag_Int package.
%(c)(1)flow chart 
% if the interpolation data are points, then use Lagrange.m function
% The Function need input the interpolation point information which is x and y,
% also it need input X, which is the x value for the interpolation
% Divided_difference to get the a=[a0,a1,...,an], which is the input for Newton_polynomial
%, then use a,x,X to get the interpolation value at X.
% if the interpolation data is a function, then use
% Lagrange_polynomial_interpolation.m .It need input lower bound a, and
% upper bound b, n many space point,function f, and target evaluation
% points X.
%(2) how I test the correctness of my overall implementation. Using
%Lagrange_polynomial_interpolation.m fuction and large enough n point to
%interpolate many different functions, and see if the interpolation line close to
%the original function. 
%% (d)
n = 5;
g = @(x) sin(x);
m = 1000;  % number of evaluation points
f= @(x) 1./(1+25*x.^4); % define the function
x = (0:2/(n-1):2)-1;
z = (0:2/(m-1):2)-1; % x location of evaluation points
zz =(0:2/(m-1):pi);
[Pf]=Lagrange_polynomial_interpolation(-1,1,f,n,z);
[Pg]=Lagrange_polynomial_interpolation(0,pi,g,n,zz);
%% Question 2
%%ploting for natural cubic spline.
i = 1;

for n = [5,11,21] % using increasing number of data points    
    x = (0:2/(n-1):4)-2; % x location of data points
    z = (0:2/(m-1):2)-1;
    y = f(x); % f(x) values
    %% evaluate S(z), the polynomial at z
    [a,b,c,d] = natural_spline_coeff(x,y);
    for j=1:length(z)
          Sz(j) =natural_spline_polynomial(a,b,c,d,x,z(j));
    end
    %% evaluate P(z), the polynomial at z
    [Pz]=Lagrange_polynomial_interpolation(-2,2,f,n,z);
    %% Ploting the Spline_polynomial, Newton_polynomial, and exact function
    figure(i); 
    plot(z,Pz);
    hold on
    plot(z,Sz);
    plot(z,f(z));
    legend('Newton','cubic','exact function','Location','northwest');

    %% compute the abs error of the interpolation at z    
    errorP(i)= norm(Pz-f(z))/sqrt(m);
    errorS(i)= norm(Sz-f(z))/sqrt(m);
    errorP1(i)=norm(Pz([m/4:3*m/4])-f(z([m/4:3*m/4])))/sqrt(m);
    i = i+1;
end
    %% plot errors
figure(4); plot([5,11,21].*log([5,11,21]),log(errorP));
figure(5); plot([5,11,21].*log([5,11,21]),log(errorS));
figure(6); plot([5,11,21].*log([5,11,21]),log(errorP1));
%% the clamped
n=21;
x = (0:2/(n-1):2)-1; 
y=f(x);
ff= @(x) -((100*x.^3)./(1+25*x.^4).^2);
yb=ff(x);
y0=ff(-1);
yn=ff(1);
for i=1:21
    yyy(3*i-2)=yb(i);
    yyy(3*i-1)=y(i);
    yyy(3*i)=yb(i);
end
z = (0:2/(m-1):2)-1;
cs = spline(x,[ yyy ]);
errorSC= norm(ppval(cs,z)-f(z))/sqrt(m);
errorSN= norm(Sz-f(z))/sqrt(m);
plot(x,y,'o',z,ppval(cs,z),'-');
%% Question 2 discussion
% I quantify the error with the infinity norm ||P(x)-f(x)||, I also tried
% to normalize it by divide sqrt(m), where m is the number of evaluation
% points. With less error, the interpolation is more accuacy.
% The effect of boundary condition on the spline interpolation.
%% functions
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

function [PX] = Newton_polynomial(b,x,X)
n = numel(x);
S = ones(size(X));
PX = b(1)*S;

for i = 2:n
        S = S.*(X-x(i-1));
        PX = PX+ b(i).*S;
end


end
function  [Lx]=Lagrange_polynomial_interpolation(a,b,f,n,j)
% a is lower bound, b is upper bound, n is the number of point that want to
% interpolate
x = (0:2/(n-1):b-a)+a; %Divide in terval [a,b] in to n-1 area, which generate n point.
y=f(x);%Get the function value
[alpha] = Divided_difference(x,y);
Lx = zeros(size(j));

[Lx] = Newton_polynomial(alpha,x,j);
end

function [a,b,c,d]=natural_spline_coeff(x,y)
n=numel(x);
a=y;
h=[];e=[];l=[];u=[];z=[];
for i=1:n-1
    h(i)=x(i+1)-x(i);
end
for i=2:n-1
    e(i)=(((a(i+1)-a(i))/h(i))-(a(i)-a(i-1))/h(i-1))*3;
end
l=1;
u=0;
z=0;
for i=2:n-1
    l(i)=2*(x(i+1)-x(i-1))-h(i-1)*u(i-1);
    u(i)=h(i)/l(i);
    z(i)=(e(i)-h(i-1)*z(i-1))/l(i);
end
l(n)=1;
z(n)=0;
c=[];b=[];d=[];
c(n)=0;
for j= n-1:-1:1
    c(j)=z(j)-u(j)*c(j+1);
    b(j)=(a(j+1)-a(j))/h(j)-h(j)*(c(j+1)+2*c(j))/3;
    d(j)=(c(j+1)-c(j))/(3*h(j));
end
end
function [Sz]=natural_spline_polynomial(a,b,c,d,x,z)
n=numel(x);
Sz=a(1);
for i=1:n-1
    if((z>x(i))&(z<=x(i+1)))
        Sz=a(i)+b(i)*(z-x(i))+c(i)*(z-x(i))^2+d(i)*(z-x(i))^3;
    end
end
end