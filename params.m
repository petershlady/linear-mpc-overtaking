% Script to store all the parameters used for finding a safe path within
% stability limits and around obstacles. Stores values in a struct. that can
% be used elsewhere.
% 
% Author: Peter Schleede

addpath('utils');

P.variable_Ux_plan  = 1;
P.driving_mode      = 1;    % 1- normal, 2 - waiting/slow
P.is_oncoming       = 1;

%% initial conditions
P.x.beta_init       = 0;
P.x.r_init          = 0;
P.x.dpsi_init       = 0;
P.x.s_init          = 0;
P.x.e_init          = .10;

%% constraints
P.con.Fmax          = 4000;
P.con.Fmin          = -P.con.Fmax;
P.con.Vmax          = 500;
P.con.Vmin          = -P.con.Vmax;

%% problem parameters
P.prob.num_states   = 5;

% time vector used at each step
P.prob.t_short      = 0.05;                 % s
P.prob.t_long       = 0.25;                 % s
P.prob.T_corr       = 10;
P.prob.T_long       = 30;
P.prob.dt           = 0.05;                 % s, actual time step we are using
P.prob.t_startup    = 1;                    % do not ever change in this 
                                            % file, used for time vector 
                                            % initialization
P.prob.num_steps    = 400;                    % number of loop iterations
P.prob.obst_mode    = 2;
P.prob.pass_side    = 'left';

%% vehicle parameters
P.veh.width         = 1.75;                 % meters
P.veh.mass          = 1200;                 % kg
P.veh.Fz            = 9.81 * P.veh.mass;    % N
P.veh.Ux_nom        = 25;
P.veh.Ux_slow       = 17.5;
P.veh.Ux            = P.veh.Ux_nom;         % m/s
P.veh.Ux_actual     = P.veh.Ux_nom;
P.veh.ax_max        = 4;
P.veh.ax_min        = -1.4;
P.veh.a             = 1.25;                 % m
P.veh.b             = 2.0;                  % m
P.veh.L             = P.veh.a + P.veh.b;    % m
P.veh.mu            = 0.75;
P.veh.Ca            = 110000;
P.veh.Ca_lin        = 80000;
P.veh.tire_mode     = 'fiala';
P.veh.Izz           = (1/12) * P.veh.mass * ((P.veh.a+P.veh.b)^2+2^2);

%% path parameters
P.path.e_buffer     = 1;                    % meters
P.path.e_max        = 6;
P.path.e_min        = -2;

%% optimization parameters
% weights for slack variables
P.opt.W_veh         = [1, 1, 1, 1]';
P.opt.W_env         = [3e2, 3e2]';

% quadratic terms
P.opt.Q_beta        = 0;
P.opt.Q_r           = 0;
P.opt.Q_dpsi        = 1000;
P.opt.Q_s           = 0;
P.opt.Q_e           = 2;
P.opt.Q             = [P.opt.Q_beta, 0, 0, 0, 0;
                       0, P.opt.Q_r, 0, 0, 0;
                       0, 0, P.opt.Q_dpsi, 0, 0;
                       0, 0, 0, P.opt.Q_s, 0;
                       0, 0, 0, 0, P.opt.Q_e];
% https://www.mathworks.com/matlabcentral/answers/...
% 324971-forming-a-block-diagonal-matrix-of-one-certain-matrix
Q_tilde             = repmat(P.opt.Q, 1, P.prob.T_long);
Q_tilde             = mat2cell(Q_tilde, size(P.opt.Q,1), repmat(size(P.opt.Q,2),1,P.prob.T_long));
P.opt.Q_tilde       = sparse(blkdiag(Q_tilde{:}));
P.opt.Q_e2          = 0;
P.opt.R             = .1;
P.opt.v_scale       = 100;

%% others
P.vis.colors        = lines(5);
P.vis.plot_trajs    = 1;
P.vis.trajectories  = 1;

%% target vehicle
target.width        = 1.75;                 % meters
target.mass         = 1200;                 % kg
target.Fz           = 9.81 * target.mass;   % N
target.Ux           = 17.5;                 % m/s
target.a            = 1.25;                 % m
target.b            = 2.0;                  % m
target.L            = target.a + target.b;  % m
target.mu           = 0.75;
target.Ca           = 80000;
target.Ca_lin       = 110000;
target.tire_mode    = 'fiala';
target.Izz          = (1/12) * target.mass * ((target.a+target.b)^2+2^2);

% initial conditions
target.beta0        = 0;
target.r0           = 0;
target.dpsi0        = 0;
target.s0           = 37.5;
target.e0           = 0;
target.x            = [target.beta0;...
                       target.r0;...
                       target.dpsi0;...
                       target.s0;...
                       target.e0];
                   
%% oncoming vehicle
oncoming.width        = 1.75;                 % meters
oncoming.mass         = 1200;                 % kg
oncoming.Fz           = 9.81 * oncoming.mass; % N
oncoming.Ux           = -17.5;                % m/s
oncoming.a            = 1.25;                 % m
oncoming.b            = 2.0;                  % m
oncoming.L            = oncoming.a + oncoming.b;  % m
oncoming.mu           = 0.75;
oncoming.Ca           = 80000;
oncoming.Ca_lin       = 110000;
oncoming.tire_mode    = 'fiala';
oncoming.Izz          = (1/12) * oncoming.mass * ((oncoming.a+oncoming.b)^2+2^2);

% initial conditions
oncoming.beta0        = 0;
oncoming.r0           = 0;
oncoming.dpsi0        = 0;
oncoming.s0           = 200;
oncoming.e0           = 4;
oncoming.x            = [oncoming.beta0;...
                       oncoming.r0;...
                       oncoming.dpsi0;...
                       oncoming.s0;...
                       oncoming.e0];
                   
%% vehicle sim parameters (for running parametric studies)

P.sim.width         = 1.75;                 % meters
P.sim.mass          = 1200;                 % kg
P.sim.Fz            = 9.81 * P.veh.mass;    % N
P.sim.Ux            = 25;                   % m/s
P.sim.a             = 1.25;                 % m
P.sim.b             = 2.0;                  % m
P.sim.L             = P.veh.a + P.veh.b;    % m
P.sim.mu            = 0.75;
P.sim.Ca            = 80000;
P.sim.Ca_lin        = 110000;
P.sim.tire_mode     = 'fiala';
P.sim.Izz           = (1/12) * P.veh.mass * ((P.veh.a+P.veh.b)^2+2^2);
