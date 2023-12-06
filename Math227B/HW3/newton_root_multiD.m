function x = newton_root_multiD(f,J,x_0,E_tol)
     x = [];
%     x_0 =[0 0]';
%     E_tol = 10E-04;
    ep = 1;
    while E_tol < ep
        g = f(0,x_0);
        Jg = J(0,x_0);
        y = x_0 - Jg^(-1)*g;
        ep = norm(y-x_0);
        x_0 = y;
    end
    x=y;
end