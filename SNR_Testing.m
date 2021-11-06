%% Environment Clearing
clc
close all
clearvars

%% Define testing variables

% Static RNG for replicability
rng(5);

% Number of sensors
n_sensors = 3;

do_random_mean = 0;

addpath('Simulation');
addpath('DatasetGen');

%% Define the means and variance parameters

% Assume the sensors observation measurements are i.i.d.

% Each sensor's default distribution without being affected by a target
% will be a standard normal distribution with sigma = 1 and mean = [1 2 3]
mean_unaffected = [1,2,3];
var_unaffected = [1 1 1];

if do_random_mean
    mean_affected = [                                        ...
        mean_unaffected(1)+rand()*random_mean_modifier(1),   ...
        mean_unaffected(2)+rand()*random_mean_modifier(2),   ...
        mean_unaffected(3)+rand()*random_mean_modifier(3)    ...
    ];%= [2 3 4];
else
    mean_affected = [2,3,4];
end

% Use a static variance and changing mean
var_affected = [1 1 1];

% Generate a vector of SNRs to test with
% Assume the SNR will be measured in dB = 10*log10(P)

% Test a SNR vector between SNR_min <= SNR <= SNR_max using n_trials test
% sequences
SNR_min = -5;
SNR_max = 10;
n_trials = 30 + 1;

% Form the SNR vector
SNR = linspace(SNR_min, SNR_max, n_trials);

means = zeros(n_trials,length(mean_unaffected));
% Determine the new means to test with depending on the desired SNR
for i = [1:n_trials]
    means(i,:) = var_unaffected .* db2pow(SNR(i)) + mean_unaffected;
end

%% Testing Scenarios

% Number of testing scenarios
n_scenarios = 100;

% Define matrices that will hold the ADDs and PFAs for the not identical
% trials
ADD_cv = zeros(n_scenarios,n_trials);
PFA_cv = zeros(n_scenarios,n_trials);

% Make a progress bar to see the simulation progress
h = waitbar(0, 'Simulation Progress');

% Perform a test scenario where the variances are not identical for the
% unaffected and affected distributions
for i = [1:n_trials]
    for j = [1:n_scenarios]
         % State sequence, observations and changepoint
        [X_B, y_B, nu_B] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
            means(i,:),var_affected);
        
        [c_ADD, c_PFA, ~, ~] = BayesianScenario(mean_unaffected, ...
            var_unaffected,  means(i,:), var_affected, y_B, nu_B, 0.995);
        
        % Check if the change was detected at all
        if isempty(c_ADD) || isempty(c_PFA)
            ADD_cv(j,i) = nan;
            PFA_cv(j,i)  = nan;
        else
            ADD_cv(j,i) = c_ADD;
            PFA_cv(j,i)  = c_PFA;
        end
        waitbar(((i*n_scenarios) + (j-1))/(n_trials*n_scenarios));
    end
end

close(h);

% % Define matrices that will hold the ADDs and PFAs for the identical
% % trials
% ADD_iv = zeros(n_scenarios,n_trials);
% PFA_iv = zeros(n_scenarios,n_trials);
% 
% 
% % Perform a test scenario where the variances for pre and post-change
% % distribution are identical with only the mean changing
% for i = [1:n_trials]
%     for j = [1:n_scenarios]
%         [ADD_iv(j,i), PFA_iv(j,i)] = BayesianScenario(mean_unaffected, ...
%             var_unaffected, mean_affected, var_affected);
%     end
% end

%% Plot the results

% Plot the PFA and log of the ADD against each other to
% determine the spread of the simulation results
figure
colororder({'b','m'})
yyaxis left
plot(SNR, mean(ADD_cv,1),'b*') % ADD
ylabel('Average Detection Delay','Interpreter','Latex')
yyaxis right
plot(SNR, mean(PFA_cv,1),'m.') % PFA
ylabel('Average Probability of False Alarm','Interpreter','Latex')
title('ADD and PFA vs. SNR')
xlabel('SNR (dB)')
set(gca, 'color', [0 0.07 0.1 0.2])
