function [x1, Ux_actual] = simulate(x, Fyf, alpha_r, K, dt, f_inv, P)
% Function to simulate one time step and update the state.
% 
% Inputs:
%   x:                  current state
%   Fyf:                front tire force
%   alpha_r:            rear tire slip angle
%   K:                  path curvature
%   dt:                 time step
%   f_inv:              f_tire inverse lookup table
%   P:                  parameter struct
% 
% Ouputs:
%   x1:                 state after time step
%   Ux_actual:          Ux after acceleration
% 
% Usage:
%   x1 = simulate(x, Fyf, alpha_r, K, dt, f_inv, P);
% 
% History:
%   Peter Schleede, 4/21/19 - Initial version with linear dynamics.
%   Peter Schleede, 4/25/19 - Updated for non-linear dynamics.
%   Peter Schleede, 5/30/19 - Added longitudinal dynamics.
%   Peter Schleede, 6/04/19 - Added variable Ux to planning

Fyr = f_tire(alpha_r, 'fiala', P);

% If we are needing to change speeds, we need to derate our acceleration
% ability by the required front tire lateral force. We assume that the
% brakes and engine only act on the front wheels.
Ux_actual = P.veh.Ux_actual;
if Ux_actual ~= P.veh.Ux
    Fx_max = sqrt((P.veh.mu*P.veh.Fz)^2 - Fyf^2);   % unsigned
    
    % find the amount of acceleration we can actually achieve
    if Ux_actual < P.veh.Ux
        ax = min(P.veh.ax_max, Fx_max/P.veh.mass);
    else
        ax = max(P.veh.ax_min, -Fx_max/P.veh.mass);
    end
    
    % update our speed
    new_Ux = Ux_actual + ax*P.prob.dt;
    if Ux_actual < P.veh.Ux && new_Ux > P.veh.Ux
        Ux_actual = P.veh.Ux;
    elseif Ux_actual > P.veh.Ux && new_Ux < P.veh.Ux
        Ux_actual = P.veh.Ux;
    else
        Ux_actual = new_Ux;
    end
    
end
P.sim.Ux = Ux_actual;

alpha_f = calc_f_tire_inv(f_inv, Fyf);
delta = x(1) + P.sim.a*x(2) / P.sim.Ux - alpha_f;
Uy = P.sim.Ux * tan(x(1));

xdot = zeros(5,1);
xdot(1) = (Fyf*cos(delta) + Fyr) / (P.sim.mass*P.sim.Ux) - x(2);
xdot(2) = (P.sim.a*Fyf*cos(delta) - P.sim.b*Fyr) / P.sim.Izz;
xdot(3) = x(2) - P.sim.Ux*K;
xdot(4) = P.sim.Ux*cos(x(3)) - Uy*sin(x(3));
xdot(5) = P.sim.Ux*sin(x(3)) + Uy*cos(x(3));

x1 = x + xdot*dt;

end