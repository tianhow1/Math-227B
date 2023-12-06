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

