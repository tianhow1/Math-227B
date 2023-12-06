clear
clc
f = @(t,x)[(16*x(2)/((1+x(2))*(1+10*x(2)))-1)*x(1);
(2-16*x(2)/((1+x(2))*(1+10*x(2))))*x(1)-x(2)];
J = @(t,x)[(16*x(2)/((1+x(2))*(1+10*x(2)))-1),((16-160*x(2)^2)/((1+x(2))*(1+10*x(2))))*x(1);
2-(16*x(2)/((1+x(2))*(1+10*x(2)))),((-16+160*x(2)^2)/((1+x(2))*(1+10*x(2))))*x(1)-1];
x = [];
x_0 =[1 1]';
E_tol = 0.000001;
x0=[1;0];
R=newton_root_multiD(f,J,x_0,E_tol);
%R is very close to 0, thus the equilibrium is (0,0), which satisfied my 
%analytical solution
x0=[1;0];
Y=runge_kutta2_multD(f,0,0.01,50,x0);
X=0:0.01:50;
plot(X,Y)