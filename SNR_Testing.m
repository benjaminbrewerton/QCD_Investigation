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
SNR_max = 5;
n_trials = 11;

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

% Create a matrix which holds the affected and unaffected variances
vars_unaffected = zeros(n_trials,n_sensors);
vars_affected = zeros(n_trials,n_sensors);
% Iterate to fill the variance matrices
for i = [1:n_trials]
    vars_unaffected(i,:) = mean_unaffected ./ db2pow(SNR(i));
    vars_affected(i,:) = mean_affected ./ db2pow(SNR(i));
end

%% Testing Scenarios

% Number of testing scenarios
n_scenarios = 5;

% Define matrices that will hold the ADDs and PFAs for the not identical
% trials
ADD_cv = zeros(n_scenarios,n_trials);
PFA_cv = zeros(n_scenarios,n_trials);

% Perform a test scenario where the variances are not identical for the
% unaffected and affected distributions
for i = [1:n_trials]
    for j = [1:n_scenarios]
        [ADD_cv(j,i), PFA_cv(j,i)] = BayesianScenario(mean_unaffected, ...
            vars_unaffected(i,:), mean_affected, vars_affected(i,:));
    end
end

% Define matrices that will hold the ADDs and PFAs for the identical
% trials
ADD_iv = zeros(n_scenarios,n_trials);
PFA_iv = zeros(n_scenarios,n_trials);


% Perform a test scenario where the variances for pre and post-change
% distribution are identical with only the mean changing
for i = [1:n_trials]
    for j = [1:n_scenarios]
        [ADD_iv(j,i), PFA_iv(j,i)] = BayesianScenario(mean_unaffected, ...
            var_unaffected, mean_affected, var_affected);
    end
end

%% Plot the results

% Plot the log of the PFA and log of the ADD against each other to
% determine the spread of the simulation results
% figure
% for i = [1:n_trials]
%    loglog(ADD_cv(:,i), PFA_cv(:,i), '*'), hold on % Changing Variance
%    %loglog(ADD_iv(:,i), PFA_iv(:,i), '.') % Identical Variance
% end
figure
hold on
for i = [1:n_trials]
   %plot(log(ADD_cv(:,i)), log(PFA_cv(:,i)), '*') % Changing Variance
   plot(log(ADD_iv(:,i)), log(PFA_iv(:,i)), '.') % Changing Variance
   %loglog(ADD_iv(:,i), PFA_iv(:,i), '.') % Identical Variance
end

