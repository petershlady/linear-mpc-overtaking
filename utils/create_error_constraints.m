function [e_max, e_min] = create_error_constraints(t, x, target, oncoming, P, mode, direction)
% Function to create the minimum and maximum allowed values for error based
% on the current state of the vehicle and target. Assumes there are at most
% two vehicles. Does not include things like width of the follower or any
% buffer.
% 
% Inputs:
%   t:                  time vector at this prediction step
%   x:                  state of the follower
%   target:             struct containing state of the target vehicle
%                       parameters
%   oncoming:           struct containing state of the oncoming vehicle
%                       parameters
%   P:                  struct containing follower parameters
%   mode:               method for representing target
%   direction:          'right'/'left', which side to pass the target on
% 
% Ouputs:
%   e_max:              vector of maximum allowed error at long time steps
%   e_min:              vector of minimum allowed error at long time steps
% 
% Notes:
%   modes:
%       1 - target is one long box that begins at its current position and
%           extends across the time horizon. very conservative
%       2 - target is considered a moving object with a fixed length that
%           is predicted across the time horizon
% 
% Usage:
%   [e_max, e_min] = create_error_constraints(t, x, target, oncoming, P, mode, direction);
% 
% History:
%   Peter Schleede, 5/13/19 - Initial version
%   Elliot Weiss, 5/30/19 - Runs across whole time horizon
%   Peter Schleede, 6/04/19 - Added variable Ux to planning

% make sure direction is valid
if ~(strcmp(direction, 'right') || strcmp(direction, 'left'))
    ME = MException('Envelopes:create_error_constraints', ...
               'direction must be ''right'' or ''left''');
    throw(ME);
end

% useful shorthand
Tc = P.prob.T_corr;

% figure out where we are and initialize error bounds
if P.variable_Ux_plan
    Ux_vect = P.veh.Ux * ones(size(t));
    if P.veh.Ux ~= P.veh.Ux_actual
        if P.veh.Ux <= P.veh.Ux_actual
            Ux_vect = max(P.veh.Ux_actual + P.veh.ax_min*t, P.veh.Ux);
        else
            Ux_vect = min(P.veh.Ux_actual + P.veh.ax_max*t, P.veh.Ux);
        end
    end
    predicted_s = t.*Ux_vect + x(4);
else
    predicted_s = t*P.veh.Ux + x(4);
end
e_max = P.path.e_max * ones(size(t));
e_min = P.path.e_min * ones(size(t));

if mode == 1
    % assume target is a long box
    target_s_min = target.x(4) - target.b;
    target_s_max = target.x(4) +...
                 target.Ux * t(end) + target.a;
    target_matches = find((predicted_s(Tc+1:end) >= target_s_min) .*...
                         (predicted_s(Tc+1:end) <= target_s_max));
elseif mode == 2
    % assume target and oncoming vehicles are moving at a fixed speed and 
    % we can predict where they will be (and they maintain some steady state error)
    
    target_predicted_s = t*target.Ux + target.x(4);
    target_s_min = target_predicted_s - target.b;
    target_s_max = target_predicted_s + target.a;
    target_matches = find((predicted_s >= target_s_min) .*...
                         (predicted_s <= target_s_max));
                     
    oncoming_predicted_s = t*oncoming.Ux + oncoming.x(4);
    oncoming_s_min = oncoming_predicted_s - oncoming.a;
    oncoming_s_max = oncoming_predicted_s + oncoming.b;
    oncoming_matches = find((predicted_s >= oncoming_s_min) .*...
                         (predicted_s <= oncoming_s_max));
  
else 
    ME = MException('Envelopes:create_error_constraints', ...
               'mode must be 1 or 2');
    throw(ME);
end

if ~isempty(target_matches)
    fprintf('WATCH OUT FOR TARGET CAR\n')
    if strcmp(direction, 'left')
        e_min(target_matches) = target.x(5) + target.width;
    else
        e_max(target_matches) = target.x(5) - target.width;
    end
end

if P.is_oncoming
    if ~isempty(oncoming_matches)
        fprintf('WATCH OUT FOR ONCOMING CAR\n')
        if strcmp(direction, 'left')
            e_max(oncoming_matches) = oncoming.x(5) - oncoming.width;
        end
    end
end

end
