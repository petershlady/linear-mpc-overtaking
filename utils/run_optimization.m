function [x, Fyf, v, sig_veh, sig_env, cvx_optval, cvx_status, alpha_r] = ...
    run_optimization(t, x_t, target, oncoming, K, prev_alpha_r, P)
% Function to perform the optimization problem. Allows the user to specify
% which vehicle speed to run at
% 
% Inputs:
%   t:                  time vector
%   x_t:                current vehicle state
%   target:             target vehicle struct
%   oncoming:           oncoming vehicle struct
%   K:                  curvature
%   prev_alpha_r:       rear tire slip angles from last time step
%   P:                  parameter struct
% 
% Ouputs:
%   x:                  trajectory vector
%   Fyf:                vector of front tire lateral forces
%   v:                  vector of slew rates
%   sig_veh:            stability slack variables
%   sig_env:            environmental slack variables
%   cvx_optval:         cost of the problem
%   cvx_status:         problem status
%   alpha_r:            rear tire slip angle
% 
% Usage:
% [x, Fyf, v, sig_veh, sig_env, cvx_optval, cvx_status, alpha_r] = ...
%     run_optimization(t, x_t, target, oncoming, K, prev_alpha_r, P);
% 
% History:
%   Peter Schleede, 5/30/19 - Initial version
%   Peter Schleede, 6/04/19 - Removed fixed time option

% set things up for stability constraints
% r_min and r_max are based on vehicle properties but beta limits are based
% on r so are introduced as constraints in optimization
alpha_r_peak    = atan2(3*P.veh.mass*9.81*P.veh.mu*P.veh.a, P.veh.Ca*P.veh.L);
r_max           = 9.81 * P.veh.mu / P.veh.Ux;
r_min           = -r_max;

% variables to help with shaping
Tc = P.prob.T_corr;
Tl = P.prob.T_long;
n_st = P.prob.num_states;
Tcm = Tc - 1;
Tlc = Tl - Tc + 1;

% get discrete matrices for this time step
[GsA, GsB1, GsB3, GlA, GlB1, GlB2, GlB3, alpha_r] = ...
    create_discrete_matrices(t, x_t, K, prev_alpha_r, P);

% check for target_vehicle
[e_max, e_min] = create_error_constraints(t, x_t, target, oncoming, P,...
                        P.prob.obst_mode, P.prob.pass_side);

%% Optimization problem

% variables:
%   x: [beta, r, dpsi, s, e]'
%       beta:   sideslip
%       r:      yaw rate
%       dpsi:   heading deviation
%       s:      position along path
%       e:      lateral deviation from path
%   F_y,f:      forward tire lateral force
%   v_k = F_yf,k - F_yf,k-1

% objective: 
%   sum over k (x_k'*Q*x_k + v_k'*R*v_k + W_veh*sig_veh + W_env*sig_env

% constraints:
%   x(k+1)      == A_k*x_k + B_k,1*F_yf,k + B_k,3 for k=1->10
%   x(k+1)      == A_k*x_k + B_k,1*F_yf,k + B_k,2*F_yf,k+1 + B_k,3 for
%                       k=11->29
%   H_veh*x_k   <= G_veh + sig_veh for all k
%   H_env*x_k   <= G_env + sig_env for k=11->30
%   F_yf,k      <= F_max for all k
%   abs(v_k)    <= v_max,k for all k

% definitions: 
%   sig_env, sig_veh are non-negative slack variables st the
%       problem is always feasible
%   W are weights st the slack variables are always kept small

% below version of cvx problem is to isolate issues
cvx_begin quiet
    variables x(5, Tl+1) v(Tl, 1) Fyf(Tl+1, 1) ...
              sig_veh(4, Tl+1) sig_env(2, Tl+1)

    % dynamics constraints
    x(:,1) == x_t;

    % short time steps
    reshape(x(:,2:Tc), P.prob.num_states*Tcm, 1) == ...
        GsA*reshape(x(:,1:Tcm), P.prob.num_states*(Tcm), 1) ...
        + GsB1*Fyf(1:Tcm) + GsB3;

    % long time steps
    reshape(x(:,Tc+1:Tl+1), P.prob.num_states*Tlc, 1) == ...
        GlA*reshape(x(:,Tc:Tl), P.prob.num_states*Tlc, 1) ...
        + GlB1*Fyf(Tc:end-1) + GlB2*Fyf(Tc+1:end) + GlB3;

    % force constraints
    Fyf <= P.con.Fmax;
    Fyf >= P.con.Fmin;

    % slew constraints
    v == Fyf(2:end) - Fyf(1:end-1);
    v <= P.con.Vmax;
    v >= P.con.Vmin;

    % stability constraints
    x(2, :) <= r_max + sig_veh(1, :);
    x(2, :) + sig_veh(2, :) >= r_min;
    x(1, :) <= alpha_r_peak + P.veh.b * x(2, :) / P.veh.Ux  + sig_veh(3, :);
    x(1, :) + sig_veh(4, :) >= -alpha_r_peak + P.veh.b * x(2, :) / P.veh.Ux;

    % environmental constraints
    x(5, :) <= e_max - (1/2)*P.veh.width - P.path.e_buffer + sig_env(1, :);
    x(5, :) + sig_env(2, :) >= e_min + (1/2)*P.veh.width + P.path.e_buffer;

    % slack variables
    sig_veh >= 0;
    sig_env >= 0;
    sig_env <= P.path.e_buffer;

    % objective
    cost = quad_form(reshape(x(:,2:end), P.prob.num_states*Tl,1), P.opt.Q_tilde)...
           + quad_form(v / P.opt.v_scale, P.opt.R)...
           + sum(P.opt.W_veh' * sig_veh)...
           + sum(P.opt.W_env' * sig_env);
    minimize(cost)
cvx_end

end