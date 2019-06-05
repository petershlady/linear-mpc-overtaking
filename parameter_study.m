% Script for running parametric study of model uncertainty for 
% overtaking project
% 
% Author: Andrew Shoats
% Date  : 05/17/2019
%
% Params: 

clear all
close all

% Tell main that parameter study is running
sim = 1;

% Load params
params;

%% Place changes to be run here

delta_m = -2000:500:2000; % vector of parameter changes

%% Run main (in loop)

for ii = 1:length(delta_m)
    
    P.sim.mass = P.veh.mass + delta_m(ii);
    
    main
    
end