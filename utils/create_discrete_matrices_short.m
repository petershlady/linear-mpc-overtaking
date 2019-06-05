function [GsA, GsB1, GsB3, alpha_r] = ...
            create_discrete_matrices_short(x, K, alpha_prev, P)
% Function to create discrete versions of the dynamics matrices at set
% times and based on a starting state. Follows Brown (2017) but only
% creates short term matrices.
% 
% Inputs:
%   t:          time vector at which to discretize
%   x:          initial state
%   K:          curvature at next time steps (assume traveling along
%               centerline at Ux and use t)
%   alpha_prev: vector of rear tire slip angles at last time step
%   P:          parameter struct
% 
% Outputs:
%   GsA:        3d matrix A for short time steps
%   GsB1:       3d matrix B1 for short time steps
%   GsB3:       3d matrix B3 for short time steps
%   alpha_r:    rear tire slip angle at this time step
% 
% Usage:
%   [GsA, GsB1, GsB3, GlA, GlB1, GlB2, GlB3, alpha_r] = ...
%             create_discrete_matrices(t, x, K, alpha_prev, P);
% 
% History:
%   Peter Schleede, 4/17/19 - Initial version
%   Peter Schleede, 4/18/19 - Added continuous and discrete, modified
%                             outputs
%   Peter Schleede, 4/20/19 - Fixed matrix exponentials and Gamma1
%   Peter Schleede, 4/21/19 - Fixed off by one error in long matrices
%   Peter Schleede, 5/08/19 - Updated dynamics to use correct rear slip
%                             angles

n_st = P.prob.num_states;

%% initialize continuous matrices
A_c = zeros(n_st, n_st, P.prob.T_long);
B_c_F = zeros(n_st, 1);
B_c_K = zeros(n_st, 1);
d_c = zeros(n_st, 1, P.prob.T_long);

% get linearized rear tire forces (based off past slip angles)
[Ca_lin, Fyr] = calculate_rear_tire_forces(alpha_prev, P);

% calculate alpha_r
alpha_r = x(1) - P.veh.b * x(2) / P.veh.Ux;

%% fill in continuous matrices
% linearized beta dot is (F_y,f+F_y,r) / (m*Ux) - r
A_c(1, 1, :) = -Ca_lin / (P.veh.mass * P.veh.Ux);
A_c(1, 2, :) = Ca_lin*P.veh.b / (P.veh.mass * P.veh.Ux^2) - 1;
B_c_F(1) = 1 / (P.veh.mass * P.veh.Ux);
d_c(1, 1, :) = (Ca_lin.*alpha_prev + Fyr) / (P.veh.mass * P.veh.Ux);

% linearized r dot is (a*F_y,f - b*F_y,r) / Izz
A_c(2, 1, :) = P.veh.b*Ca_lin / P.veh.Izz;
A_c(2, 2, :) = -P.veh.b^2*Ca_lin / (P.veh.Ux * P.veh.Izz);
B_c_F(2) = P.veh.a / P.veh.Izz;
d_c(2, 1, :) = -(P.veh.b*Fyr + Ca_lin.*alpha_prev) / P.veh.Izz;

% linearized d_psi dot is r - Ux*K(s)
A_c(3, :, :) = repmat([0, 1, 0, 0, 0], 1, 1, P.prob.T_long);
B_c_K(3) = -P.veh.Ux;

% linearized s dot is Ux
d_c(4, 1, :) = P.veh.Ux;

% linearized e dot is Ux * d_psi + Ux * beta
A_c(5, :, :) = repmat([P.veh.Ux, 0, P.veh.Ux, 0, 0], 1, 1, P.prob.T_long);
%% discretize matrices
% Tc = P.prob.T_corr;
Tl = P.prob.T_long;

% pre-allocate
G_short_disc = zeros(n_st+2, n_st+2, Tl);

% fill in matrices
% short term
for i=1:Tl
    G_short = [A_c(:,:,i), B_c_F, d_c(:,:,i) + B_c_K*K(i);...
                 zeros(2,n_st+2)];
    G_short_disc(:,:,i) = expm(G_short * P.prob.t_short);
end

% extract state update matrices
GsA = G_short_disc(1:n_st, 1:n_st, :);
GsB1 = G_short_disc(1:n_st, n_st+1, :);
GsB3 = G_short_disc(1:n_st, n_st+2, :);

end