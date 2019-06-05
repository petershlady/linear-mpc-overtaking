% Script to perform the optimization described in "Safe Driving Envelopes
% for Path Tracking in Autonomous Vehicles".
% 
% Author: Peter Schleede
clear

% get parameter struct
run('params');

% initial condition
x0 = [P.x.beta_init, P.x.r_init, P.x.dpsi_init, P.x.s_init, P.x.e_init]';

%% Initialize variables for the problem
% previous alpha_rear and delta will be updated but for now set to zeros
prev_alpha_r    = zeros(P.prob.T_long, 1);

% set things up for stability constraints
% r_min and r_max are based on vehicle properties but beta limits are based
% on r so are introduced as constraints in optimization
alpha_r_peak    = atan2(3*P.veh.mass*9.81*P.veh.mu*P.veh.a, P.veh.Ca*P.veh.L);
r_max           = 9.81 * P.veh.mu / P.veh.Ux;
r_min           = -r_max;

% create f_tire inverse lookup table
f_inv           = create_f_tire_inv(P);

% set fixed curvature for now
K               = 0*ones(P.prob.T_long+1, 1);
% K               = 0.005*ones(P.prob.T_long+1, 1);
% K              = -0.005*ones(P.prob.T_long+1, 1);

% pre-allocate spaces for data
t_hist          = zeros(1, P.prob.num_steps+1);
x_hist          = zeros(P.prob.num_states, P.prob.num_steps+1);
target_hist     = zeros(P.prob.num_states, P.prob.num_steps+1);
oncoming_hist   = zeros(P.prob.num_states, P.prob.num_steps+1);
F_hist          = zeros(1, P.prob.num_steps);
cost_hist       = zeros(5, P.prob.num_steps);
traj_hist       = zeros(P.prob.num_states, P.prob.T_long, P.prob.num_steps);
sigenv1_hist    = zeros(P.prob.T_long+1, P.prob.num_steps);
sigenv2_hist    = zeros(P.prob.T_long+1, P.prob.num_steps);
Ux_hist         = zeros(1, P.prob.num_steps+1);

% set x as x0
x_t             = x0;
x_hist(:,1)     = x_t;
target_hist(:,1) = target.x;
oncoming_hist(:,1) = oncoming.x;
Ux_hist(1)      = P.veh.Ux_actual;

%% Main Loop
%   At each timestep, need to:
%       1. Create time vector
%       2. Create discrete linear system
%       3. Run optimization problem
%       4. Update variables as necessary
tic
for k=1:P.prob.num_steps
    % create time vector
    t = create_time_vector(P);
    if P.prob.t_startup
        P.prob.t_startup = 0;
    end
    
    if P.driving_mode == 1
        P.veh.Ux = P.veh.Ux_nom;
        [x, Fyf, v, sig_veh, sig_env, cvx_optval, cvx_status, alpha_r] = ...
            run_optimization(t, x_t, target, oncoming, K, prev_alpha_r, P);
    end
    
    % when we cannot safely pass, we want to move to mode 2
    if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Inaccurate/Solved')
        P.driving_mode = 2;
    end
    
    if P.driving_mode == 2
        P.veh.Ux = P.veh.Ux_slow;
        [x, Fyf, v, sig_veh, sig_env, cvx_optval, cvx_status, alpha_r] = ...
        run_optimization(t, x_t, target, oncoming, K, prev_alpha_r, P);
    end
    
    % if the oncoming vehicle is past, return to original problem
    if (x_t(4) - oncoming.x(4)) > 10
        P.driving_mode = 1;
    end
    
    % simulate one time step
    [x_t, P.veh.Ux_actual] = simulate(x_t, Fyf(1), alpha_r, K(1), P.prob.dt, f_inv, P);
    target.x(4) = target.x(4) + P.prob.dt*target.Ux;
    oncoming.x(4) = oncoming.x(4) + P.prob.dt*oncoming.Ux;
    
    % display vehicle positions
    fprintf('Our position: %f\n', x_t(4))
    fprintf('Target position: %f\n', target.x(4))
    if P.is_oncoming
        fprintf('Oncoming position: %f\n', oncoming.x(4))
    end
    
    % update previous alpha_r vector
    prev_alpha_r = (x(1,2:end) - P.veh.b*x(2,2:end) / P.veh.Ux)';
    
    % save observed results
    t_hist(k+1)             = k*P.prob.t_short;
    x_hist(:, k+1)          = x_t;
    target_hist(:, k+1)     = target.x;
    oncoming_hist(:, k+1)   = oncoming.x;
    F_hist(:, k)            = Fyf(1);
    cost_hist(1, k)         = cvx_optval;
    cost_hist(2, k)         = sum(x(3,2:end).^2*P.opt.Q_dpsi);
    cost_hist(3, k)         = sum(x(5,2:end).^2*P.opt.Q_e);
    cost_hist(4, k)         = sum((v/P.opt.v_scale).^2*P.opt.R);
    cost_hist(5, k)         = sum(P.opt.W_veh'*sig_veh);
    traj_hist(:, :, k)      = x(:, 2:end);
    sigenv1_hist(:,k)       = sig_env(1,:);
    sigenv2_hist(:,k)       = sig_env(2,:);
    Ux_hist(k+1)            = P.veh.Ux_actual;
    
    % let user know about progress
    fprintf('Finished step %d/%d, %.3f seconds of simulation time\n',...
                k, P.prob.num_steps, k*P.prob.dt);
    fprintf('Current speed: %f\n', P.veh.Ux_actual)
    fprintf('CVX status: %s\n', cvx_status);
    
end
toc
%% results visualization
title_font = 18;
label_font = 16;

figure(1), clf
hold on
% trajectories
if P.vis.plot_trajs
    for i=1:P.prob.num_steps / P.vis.trajectories
        traj_i = traj_hist(:, :, i*P.vis.trajectories);
        plot(-traj_i(5,:), traj_i(4,:), 'Color', P.vis.colors(3,:))
    end
end
% s and e
plot(zeros(1,P.prob.num_steps+1), x_hist(4,:), 'Color', P.vis.colors(2,:),...
        'LineStyle', '--')
plot(-x_hist(5, :), x_hist(4, :), 'Color', P.vis.colors(1,:), 'LineWidth', 1.5)
xlim([-P.path.e_max,-P.path.e_min])
xlabel('e (m)', 'FontSize', label_font)
ylabel('s (m)', 'FontSize', label_font)
title('Error along path', 'FontSize', title_font)

% heading error
figure(2), clf
plot(t_hist, x_hist(3, :), 'Color', P.vis.colors(1,:))
xlabel('time', 'FontSize', label_font)
ylabel('\Delta\Psi', 'FontSize', label_font)
title('Heading error', 'FontSize', title_font)

% objective
figure(3), clf
plot(t_hist(2:end), cost_hist(1,:), 'Color', P.vis.colors(1,:))
hold on
plot(t_hist(2:end), cost_hist(2,:), 'Color', P.vis.colors(2,:))
plot(t_hist(2:end), cost_hist(3,:), 'Color', P.vis.colors(3,:))
plot(t_hist(2:end), cost_hist(4,:), 'Color', P.vis.colors(4,:))
xlabel('time', 'FontSize', label_font)
ylabel('Cost', 'FontSize', label_font)
title('Costs', 'FontSize', title_font)
legend('Objective cost', 'Heading error cost', 'Lateral error cost', 'Slew rate cost')

% steering angles
delta = rad2deg(x_hist(1, 1:end-1) + P.veh.a*x_hist(2, 1:end-1) /...
            P.veh.Ux - calc_f_tire_inv(f_inv, F_hist));
figure(4), clf
plot(t_hist(1:end-1), delta, 'Color', P.vis.colors(1,:))
xlabel('time', 'FontSize', label_font)
ylabel('\delta', 'FontSize', label_font)
title('Steering angles in degrees', 'FontSize', title_font)

% check stability
beta_max_u = alpha_r_peak + P.veh.b*r_max / P.veh.Ux;
beta_max_l = alpha_r_peak + P.veh.b*r_min / P.veh.Ux;
beta_min_l = -(alpha_r_peak + P.veh.b*r_max / P.veh.Ux);
beta_min_u = -(alpha_r_peak + P.veh.b*r_min / P.veh.Ux);
r = r_min:.01:r_max;
beta = alpha_r_peak + P.veh.b*r ./ P.veh.Ux;
beta_min = -alpha_r_peak +P.veh.b*r ./ P.veh.Ux;

figure(5), clf
hold on
grid on
plot([beta_min_l,beta_max_l], [r_min,r_min], 'b')
plot([beta_min_u,beta_max_u], [r_max,r_max], 'b')
plot(beta, r, 'b')
plot(beta_min, r, 'b')
scatter(x_hist(1, :), x_hist(2, :), 'r', 'x')
xlabel('\beta', 'FontSize', label_font)
ylabel('r', 'FontSize', label_font)
title('Stability envelope', 'FontSize', title_font)

figure(6), clf
plot(t_hist, Ux_hist)
title('Actual vehicle speed', 'FontSize', title_font)
xlabel('Time (s)', 'FontSize', label_font)
ylabel('U_x (m/s)', 'FontSize', label_font)
