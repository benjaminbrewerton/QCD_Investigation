%% Cleanup

clc
close all
clear

% Add required paths
addpath('Simulation')
addpath('DatasetGen')

%% Testing Variables

% Script for testing the MTFA

% Change the simulation command depending on what MTFA is required.

% Define the number of trials to conduct
n_trials = 200;

% Define a mean and variance to test from
var_mu = [ 1 1 1 ]; var_nu = [ 1 1 1 ];
mean_mu = [ 1 2 3 ]; mean_nu = [ 2 3 4 ];

%% Execution for varying mean
% 
% % Test a SNR vector between SNR_min <= SNR <= SNR_max using n_trials test
% % sequences
% SNR_min = -5;
% SNR_max = 5;
% n_snr = 10 + 1;
% 
% % Form the SNR vector
% SNR = linspace(SNR_min, SNR_max, n_snr);
% 
% means = zeros(n_snr,length(mean_mu));
% % Determine the new means to test with depending on the desired SNR
% for i = [1:n_snr]
%     means(i,:) = var_mu .* db2pow(SNR(i)) + mean_mu;
% end
% 
% % Make a progress bar to see the simulation progress
% u = waitbar(0, 'Simulation Progress');
% 
% % Iterate and take the maximum value of the test statistic
% stat = zeros(1,n_snr);
% for j = [1:n_snr]
%     for i = [1:n_trials]
%        % Get a scenario for the simulation
%        [~,y,lambda] = simulate_random_scenario(mean_mu, var_mu, means(j,:), ...
%            var_nu, 0);
% 
%        % Define a threshold
%        [~,~,~,stat_cur] = BayesianScenario(mean_mu, var_mu, means(j,:), ...
%            var_nu, y, lambda);
% 
%        % Compare the stat with known maximum
%        stat_max = max(stat_cur(2,:));
%        %stat_max = max(stat_cur);
%        stat(j) = max(stat(j), stat_max);
% 
%        % Update the progress bar
%        waitbar((((j-1)*n_trials) + i)/(n_snr*n_trials));
%     end
% end

%% Static Mean

lambda_sum = 0;

stat = 0;
for i = [1:n_trials]
   % Get a scenario for the simulation
   [X,y,lambda] = simulate_random_scenario(mean_mu, var_mu, mean_nu, ...
       var_nu, 0);
   
   lambda_sum = lambda_sum + lambda;
   
   % Define a threshold
   [~,~,~,stat_cur] = BayesianScenario(mean_mu, var_mu, mean_nu, ...
       var_nu, y, lambda);

   % Compare the stat with known maximum
   stat_max = max(stat_cur(2,:));
   %stat_max = max(stat_cur);
   stat = max(stat, stat_max);

   % Update the progress bar
   waitbar(i/n_trials);
end

lambda_avg = lambda_sum / n_trials;