%% Environment Clearing
clc
close all
clearvars

%% Define testing variables

% Static RNG for replicability
rng(5);

% Number of sensors
n_sensors = 3;

%% Produce a set of statistical SNRs for testing

% Assume the SNR will be measured in dB = 10*log10(P)

% Test a SNR vector between SNR_min <= SNR <= SNR_max using n_trials test
% sequences
SNR_min = -5;
SNR_max = 15;
n_trials = 61;

% Form the SNR vector
SNR = linspace(SNR_min, SNR_max, n_trials);

%% Produce Means and Variances from the SNR vector

% Assume static means for the affected and unaffected distribution
% sequences. The variance will be changed according to the SNR
mean_unaffected = [1 2 3];
mean_affected = [2 3 4];

% Create some static variances as well for testing
var_unaffected = [1 1 1];
var_affected = [1 1 1];

% Define the limits of Variance
% Unaffected
a = min(var_unaffected); % Top limit
b = max(var_unaffected); % Bottom limit

% Affected
c = min(var_affected); % Top limit
d = max(var_affected); % Bottom limit

% Create a matrix which holds the affected and unaffected variances
vars_unaffected = zeros(n_trials,n_sensors);
vars_affected = zeros(n_trials,n_sensors);
% Iterate to fill the variance matrices
for i = [1:n_trials]
    vars_unaffected(i,:) = mean_unaffected ./ db2pow(SNR(i));
    vars_affected(i,:) = mean_affected ./ db2pow(SNR(i));
    % Use a random variance to simulate with
%     vars_unaffected(i,:) = a + (b-a) .* rand(n_sensors,1);
%     vars_affected(i,:) = c + (d-c) .* rand(n_sensors,1);
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
        [c_ADD, c_PFA, ~, ~] = BayesianScenario(vars_unaffected(i,:), ...
            var_unaffected,  vars_affected(i,:), var_affected);
        
        % Check if the change was detected at all
        if isempty(c_ADD) || isempty(c_PFA)
            ADD_cv(j,i) = nan;
            PFA_cv(j,i)  = nan;
        else
            ADD_cv(j,i) = c_ADD;
            PFA_cv(j,i)  = c_PFA;
        end
        waitbar(((i*n_scenarios) + j)/(n_trials*n_scenarios));
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
title('ADD and PFA vs. distribution SNR')
xlabel('SNR (dB)')
set(gca, 'color', [0 0.07 0.1 0.2])
