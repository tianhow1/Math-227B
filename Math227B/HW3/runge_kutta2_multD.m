function Y=runge_kutta2_multD(F,t0,h,tf,y0)
y=y0;
Y=y;
for t=t0:h:tf-h
         k1 = F(t,y);
         k2 = F(t+h, y+h*k1);
         y = y + h*(k1 + k2)/2;
         Y=[Y,y];
end
end
