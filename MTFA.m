%% Cleanup

clc
close all
clear

%% Execution

% Script for testing the MTFA

% Change the simulation command depending on what MTFA is required.

% Define the number of trials to conduct
n_trials = 2500;

% Define a mean and variance to test from
var_mu = [ 1 1 1 ]; var_nu = [ 1 1 1 ];
mean_mu = [ 1 2 3 ]; mean_nu = [ 2 3 4 ];

% Make a progress bar to see the simulation progress
u = waitbar(0, 'Simulation Progress');

% Iterate and take the maximum value of the test statistic
stat = 0;
for i = [1:n_trials]
   % Get a scenario for the simulation
   [~,y,lambda] = simulate_ISD_scenario(mean_mu, var_mu, mean_nu, ...
       var_nu, 0);
   
   % Define a threshold
   [~,~,~,stat_cur] = ISD_Random(mean_mu, var_mu, mean_nu, ...
       var_nu, y, lambda);
   
   % Compare the stat with known maximum
   stat_max = max(stat_cur(2,:));
   %stat_max = max(stat_cur);
   stat = max(stat, stat_max);
   
   % Update the progress bar
   waitbar(i/n_trials);
end