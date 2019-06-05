function xdot = xfun(t, x, Fyf, K, P)

alpha_r = x(1) - P.veh.b*x(2)/P.veh.Ux;
Fyr = f_tire(alpha_r, 'linear', P);

A = [0, -1, 0,0,0;...
     0,0,0,0,0;...
     0,1,0,0,0;...
     0,0,0,0,0;...
     P.veh.Ux,0,P.veh.Ux,0,0];
b1 = [1/(P.veh.mass*P.veh.Ux), P.veh.a/P.veh.Izz,0,0,0]';
b2 = [0,0,-P.veh.Ux,0,0]';
b3 = [Fyr/(P.veh.mass*P.veh.Ux), -P.veh.b*Fyr/P.veh.Izz,0,20,0]';

big_guy = [A, b1, b2, b3, zeros(5, 3);...
    zeros(3,8), eye(3)/.2;...
    zeros(3,11)];
out.big_guy = big_guy;
out.disc = expm(big_guy*.2);
save('big_guy.mat', 'out')

xdot = A*x+b1*Fyf+b2*K+b3;
end