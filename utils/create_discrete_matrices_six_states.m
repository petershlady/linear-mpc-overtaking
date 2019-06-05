function [GsA, GsB1, GsB3, GlA, GlB1, GlB2, GlB3, alpha_r] = ...
            create_discrete_matrices_six_states(t, x, K, alpha_prev, Fyf_prev, x_pred, P)
% Function to create discrete versions of the dynamics matrices at set
% times and based on a starting state with variable speed.
% 
% Inputs:
%   t:          time vector at which to discretize
%   x:          initial state
%   K:          curvature at next time steps (assume traveling along
%               centerline at Ux and use t)
%   alpha_prev: vector of rear tire slip angles at last time step
%   Fyf_prev:   front tire lateral force calculated across previous time
%               horizon
%   P:          parameter struct
% 
% Outputs:
%   GsA:        block matrix A for short time steps
%   GsB1:       block matrix B1 for short time steps
%   GsB3:       block matrix B3 for short time steps
%   GlA:        block matrix A for long time steps
%   GlB1:       block matrix B1 for long time steps
%   GlB2:       block matrix B2 for long time steps
%   GlB3:       block matrix B3 for long time steps
%   alpha_r:    rear tire slip angle at this time step
% 
% Usage:
%   [GsA, GsB1, GsB3, GlA, GlB1, GlB2, GlB3, alpha_r] = ...
%             create_discrete_matrices(t, x, K, alpha_prev, P);
% 
% History:
%   Peter Schleede, 5/17/19 - Initial version

n_st = P.prob.num_states;

%% initialize continuous matrices
A_c = zeros(n_st, n_st, P.prob.T_long);
B_c_F = zeros(n_st, P.prob.T_long);
B_c_K = zeros(n_st, 1);
d_c = zeros(n_st, P.prob.T_long);

% get linearized rear tire forces (based off past slip angles)
[Ca_lin, Fyr] = calculate_rear_tire_forces(alpha_prev, P);

% calculate alpha_r
alpha_r = x(1) - P.veh.b * x(2) / P.veh.Ux;

%% fill in continuous matrices
% constant terms from linearization
betadot0 = (Ca_lin .* alpha_prev + Fyr) / (P.veh.mass * x_pred(6,:));
dUx0 = (1 ./ (P.veh.mass * x_pred(6,:).^2)) .* ...
       (Ca_lin.*(2*P.veh.b*x_pred(2,:)./x_pred(6,:)-x_pred(1,:)) + ...
       Fyf_prev' +Ca_lin.*alpha_prev+Fyr);
A_c(1, 1, :) = -Ca_lin / (P.veh.mass * x_pred(6,:));
A_c(1, 2, :) = Ca_lin*P.veh.b / (P.veh.mass * x_pred(6,:).^2) - 1;
A_c(1, 6, :) = -(1 ./ (P.veh.mass * x_pred(6,:).^2)) .* ...
       (Ca_lin.*(2*P.veh.b*x_pred(2,:)./x_pred(6,:)-x_pred(1,:)) + ...
       Fyf_prev' +Ca_lin.*alpha_prev+Fyr);
B_c_F(1, :) = 1 ./ (P.veh.mass * x_pred(6,:));
d_c(1, :) = betadot0 + dUx0;

% linearized r dot is (a*F_y,f - b*F_y,r) / Izz
rdot0 = -(P.veh.b*Fyr + Ca_lin.*alpha_prev) / P.veh.Izz;
dUx0 = -(P.veh.b^2 * Ca_lin .* x_pred(2,:) .* x_pred(6,:)) ./ ...
       (P.veh.Izz * x_pred(6,:).^2);
A_c(2, 1, :) = P.veh.b*Ca_lin / P.veh.Izz;
A_c(2, 2, :) = -P.veh.b^2*Ca_lin / (x_pred(6,:) * P.veh.Izz);
B_c_F(2, :) = P.veh.a / P.veh.Izz;
d_c(2, :) = rdot0 + dUx0;

% linearized d_psi dot is r - Ux*K(s)
%% TODO this one for Ux
A_c(3, :, :) = repmat([0, 1, 0, 0, 0], 1, 1, P.prob.T_long);
B_c_K(3) = -P.veh.Ux;

% linearized s dot is Ux
%% TODO switch to A matrix
d_c(4, :) = P.veh.Ux;

% linearized e dot is Ux * d_psi + Ux * beta
%% TODO switch to A matrix
A_c(5, :, :) = repmat([P.veh.Ux, 0, P.veh.Ux, 0, 0], 1, 1, P.prob.T_long);

%% discretize matrices
Tc = P.prob.T_corr;
Tl = P.prob.T_long;
n_st = P.prob.num_states;
Tcm = Tc - 1;
Tlc = Tl - Tc + 1;

% pre-allocate
G_short_disc = zeros(n_st+2, n_st+2, Tc-1);
G_long_disc = zeros(n_st+6, n_st+6, Tl-Tc+1);

% fill in matrices
% short term
for i=1:Tc-1
    G_short = [A_c(:,:,i), B_c_F(:, i), d_c(:,i) + B_c_K*K(i);...
                 zeros(2,n_st+2)];
    G_short_disc(:,:,i) = expm(G_short * P.prob.t_short);
end

% correlation (use first order hold and save in long term)
dtcorr = t(Tc+1) - t(Tc);
G_corr = [A_c(:,:,Tc), B_c_F(:, i), B_c_K, d_c(:,Tc), zeros(n_st,3);...
                zeros(3,n_st+3), (1/dtcorr)*eye(3);
                zeros(3,n_st+6)];
G_long_disc(:,:,1) = expm(G_corr * dtcorr);

% long term
% first one has a dt relative to t_corr
dtlong_first = t(Tc+2) - t(Tc+1);
G_long = [A_c(:,:,Tc+1), B_c_F(:, i), B_c_K, d_c(:,Tc+1), zeros(n_st,3);...
                zeros(3,n_st+3), (1/dtlong_first)*eye(3);
                zeros(3,n_st+6)];
G_long_disc(:,:,2) = expm(G_long * dtlong_first);

for i=3:(Tl-Tc)+1
    G_long = [A_c(:,:,i+Tc-1), B_c_F(:, i), B_c_K, d_c(:,i+Tc-1), zeros(n_st,3);...
                zeros(3,n_st+3), (1/P.prob.t_long)*eye(3);
                zeros(3,n_st+6)];
    G_long_disc(:,:,i) = expm(G_long * P.prob.t_long);
end

% extract state update matrices
GsA = G_short_disc(1:n_st, 1:n_st, :);
GsB1 = G_short_disc(1:n_st, n_st+1, :);
GsB3 = G_short_disc(1:n_st, n_st+2, :);

% convert into block diagonal form
GsA_cell = mat2cell(reshape(GsA, n_st, n_st*Tcm), n_st, repmat(n_st, 1, Tcm));
GsA = blkdiag(GsA_cell{:});
GsB1_cell = mat2cell(reshape(GsB1, n_st, Tcm), n_st, repmat(1, 1, Tcm));
GsB1 = blkdiag(GsB1_cell{:});
GsB3 = reshape(GsB3, n_st*Tcm, 1);

G_long_line = zeros(n_st, n_st+3, Tl-Tc+1);
for i=1:Tl-Tc+1
    Gamma1 = G_long_disc(1:n_st, n_st+1:n_st+3, i);
    Gamma2 = G_long_disc(1:n_st, n_st+4:n_st+6, i);
    G_long_line(:, 1:n_st, i) = G_long_disc(1:n_st, 1:n_st, i);
    G_long_line(:, n_st+1, i) = Gamma1(:, 1) - Gamma2(:, 1);
    G_long_line(:, n_st+2, i) = Gamma2(:, 1);
    G_long_line(:, n_st+3, i) = (Gamma1(:, 2:3) - Gamma2(:, 2:3)) * [K(i+Tc-1); 1] + ...
                                Gamma2(:, 2:3) * [K(i+Tc); 1];
end
GlA = G_long_line(1:n_st, 1:n_st, :);
GlB1 = G_long_line(1:n_st, n_st+1, :);
GlB2 = G_long_line(1:n_st, n_st+2, :);
GlB3 = G_long_line(1:n_st, n_st+3, :);

% convert into block diagonal form
GlA_cell = mat2cell(reshape(GlA, n_st, n_st*Tlc), n_st, repmat(n_st, 1, Tlc));
GlA = blkdiag(GlA_cell{:});
GlB1_cell = mat2cell(reshape(GlB1, n_st, Tlc), n_st, repmat(1, 1, Tlc));
GlB1 = blkdiag(GlB1_cell{:});
GlB2_cell = mat2cell(reshape(GlB2, n_st, Tlc), n_st, repmat(1, 1, Tlc));
GlB2 = blkdiag(GlB2_cell{:});
GlB3 = reshape(GlB3, n_st*Tlc, 1);

end