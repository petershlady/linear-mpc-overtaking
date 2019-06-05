% Script to animate target and overtaking vehicle, along with overtaking
% vehicle trajectory and lateral error as a function of time
%
% Author: Andrew Shoats
% Date  : 05/12/2019
%
% Changes to be made: 
% - Know time as a function of index, so we can know location as function
% of time (Fixed 05/13/2019)
% - Make changes for moving obstacles (Fixed 05/14/2019)
% - Verify front/rear axle animations done correctly
% - Fix "Mirror image" model in top plot (Fixed 05/13/2019)
% - Fix dynamic nature of line clear. (Line 126~) (Believe this is fixed)
%
% Additional to-do:
% - Add tires (Added 05/14/2019)
% - Figure out best way to visualize moving object (Added 05/14/2019)
% - Make animation scrollable or allow frame-by-frame
% - Add oval plot of track
% - Turn vehicle creation into function?

% Necessary vehicle params
a = P.veh.a;
b = P.veh.b;
L = P.veh.L;
w = P.veh.width;
d = 0.6858; % diameter of tires, based on common 27 in tire and wheel package

% Choose plotting mode (1 == vehicle-centered, 2 == entire-path)
mode = 2;

% Choose recording mode (1 == make GIF, else do not make GIF)
recording = 0;

% Choose if optimal path is shown (1 == show opt path, else do not)
show_opt_path = 1;


figure();

% Optimal Path
opt_path_e = reshape(traj_hist(5,1,:),1,length(traj_hist(5,1,:)));
opt_path_s = reshape(traj_hist(4,1,:),1,length(traj_hist(4,1,:)));

for ii = 1:1:length(traj_hist(1,1,:))
    
    % Overtaking Vehicle
    
    dpsi    = x_hist(3,ii); % heading deviation from path
    delta_i = deg2rad(delta(ii));    % steering input angle
    
    x_com = x_hist(4,ii); % s-location (distance along path) of center-of-mass
    x_fa  = x_com + a*cos(dpsi);    % s-location of front-axle center
    x_ra  = x_com - b*cos(dpsi);    % s-location of rear-axle center
    
    y_com = x_hist(5,ii);        % e-location of center-of-mass
    y_fa  = y_com + a*sin(dpsi); % e-location of front-axle center
    y_ra  = y_com - b*sin(dpsi); % e-location of rear-axle center

    traj_s_com = traj_hist(4,:,ii);
    traj_e_com = traj_hist(5,:,ii);
    
    % Target Vehicle
    
    dpsi_t  = target_hist(3,ii);
    delta_t = 0; % Zero steering input to target vehicle, will need to be accounted for in future iterations
    
    xt_com  = target_hist(4,ii);
    xt_fa   = xt_com + a*cos(dpsi_t); % Assumes vehicles are same size
    xt_ra   = xt_com - b*cos(dpsi_t);
    
    yt_com  = target_hist(5,ii);
    yt_fa   = yt_com + a*sin(dpsi_t);
    yt_ra   = yt_com - b*sin(dpsi_t);
    

    h1 = subplot(2,1,1);
    
        % Plot COM of overtaking-vehicle
        p1_1 = plot(y_com, x_com,'ro');
        hold on
    
        % Plot  center-line of overtaking-vehicle
        p1_2 = plot([y_ra,y_fa],[x_ra,x_fa],'b','LineWidth',1);
        
        % Plot tires of overtaking-vehicle
        p1_3 = plot([y_fa - d/2*sin(delta_i+dpsi), y_fa + d/2*sin(delta_i+dpsi)],...
                    [x_fa + d/2*cos(delta_i+dpsi), x_fa - d/2*cos(delta_i+dpsi)],...
                    'k','LineWidth',2); % Front tire
        p1_4 = plot([y_ra + d/2*sin(dpsi), y_ra - d/2*sin(dpsi)],...
                    [x_ra + d/2*cos(dpsi), x_ra - d/2*cos(dpsi)],...
                    'k','LineWidth',2); % Rear tire
    
        if show_opt_path == 1
            
            % Plot optimal path
            p1_5 = plot(opt_path_e, opt_path_s,'--b');
            
        else
            % Plot MPC trajectory plan at each time-step
            p1_5 = plot(traj_e_com, traj_s_com,'kx');
            
        end
       
    
        % Plot obstacle (target-vehicle) position (Need to fix)
        p1_6 = plot(yt_com, xt_com, 'bo');
        
        p1_7 = plot([yt_ra,yt_fa],[xt_ra,xt_fa],'k','LineWidth',1);
        
        p1_8 = plot([yt_fa - d/2*sin(delta_t+dpsi_t), yt_fa + d/2*sin(delta_t+dpsi_t)],...
                    [xt_fa + d/2*cos(delta_t+dpsi_t), xt_fa - d/2*cos(delta_t+dpsi_t)],...
                    'k','LineWidth',2); % Front tire
        p1_9 = plot([yt_ra + d/2*sin(dpsi_t), yt_ra - d/2*sin(dpsi_t)],...
                    [xt_ra + d/2*cos(dpsi_t), xt_ra - d/2*cos(dpsi_t)],...
                    'k','LineWidth',2); % Rear tire
        
    % Plot road center line (Maybe fix this so the y-lim is exact value of 
    % end of road).
        p1_10 = plot([0,0],[0,max(400*pi)],'--m');
    
    % Plot params
        xlim([-10,10])
        xlabel('e, Lateral Error (m)')
        ylabel('s, Position on Path (m)')
    
        if mode == 1 % Plot vehicle-centered
        
            ylim([x_com - 5 , x_com + 15])
        
        elseif mode == 2 % Plot 
        
            ylim([0, 75*pi])

        end
    
        if show_opt_path == 1
            legend([p1_1,p1_6,p1_5],'COM Overtaking Vehicle','COM Lead Vehicle', 'Optimal Trajectory (MPC)') 
            hold off
            
        else
            legend([p1_1,p1_6,p1_5],'COM Overtaking Vehicle','COM Lead Vehicle', 'Planned Trajectory (MPC)') 
            hold off
            
        end
         
        
            % Make animation into GIF
    % (https://www.mathworks.com/matlabcentral/fileexchange/63239-gif)
    % Credit to: Chad A. Greene UTIG, June 2017
    
    if recording == 1 % Check if recording option is selected
        
        if ii == 1
        
            gif('overtaking_animation.gif')
        
        else
        
            gif
        
        end
        
    end
        
    h2 = subplot(2,1,2);
    
        p2_1 = plot((ii-1)*P.prob.dt,y_com,'rx');
        hold on
        p2_2 = plot((ii-1)*P.prob.dt, traj_e_com(1),'kx');
        % Plot road center line
        p2_3 = plot([0,length(traj_hist(1,1,:))-1],[0,0],'--m');
        % Need changes from moving obstacle
%         rectangle('Position',[4*(obstacle_sr)/(P.prob.dt*P.veh.Ux), obstacle_el, ...
%              4*(obstacle_sf - obstacle_sr)/(P.veh.Ux), (obstacle_er - obstacle_el)],'HandleVisibility','off');
%     
%         plot(4*(obstacle_sf+obstacle_sr)/(P.veh.Ux*2),(obstacle_el+obstacle_er)/2,'bo');
        xlim([1,(length(traj_hist(1,1,:)-1)*P.prob.dt)])
        ylim([min(x_hist(5,:))-0.5,max(x_hist(5,:))+0.5])
        % Fix this legend when target vehicle is included
        legend([p2_1,p2_2],'Actual Lateral Error','Planned Lateral Error') %,'Target Vehicle') 
        xlabel('Time (s)')
        ylabel('Lateral Error (m)')
        
    pause(0.00001)
    
    if ii ~= length(traj_hist(1,1,:)) 
        
        cla(h1);
        
    end
    
end
    