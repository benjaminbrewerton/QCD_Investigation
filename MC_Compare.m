% This script is used to test the performance of the Bayesian and CUSUM
% techniques with a specific set of means and variances iaw a SNR bound
close all
clear
clc

% SET ENVIRONMENT SETTINGS
do_random_mean = 0;
random_mean_modifier = [ 1 1 1 ];

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
SNR_max = 5;
n_trials = 10 + 1;

% Form the SNR vector
SNR = linspace(SNR_min, SNR_max, n_trials);

means = zeros(n_trials,length(mean_unaffected));
% Determine the new means to test with depending on the desired SNR
for i = [1:n_trials]
    means(i,:) = var_unaffected .* db2pow(SNR(i)) + mean_unaffected;
end

%% Begin iterating around each SNR and running a number of trials
n_scenarios = 150;

% Store the results in matrices
results_CUSUM_B = zeros(n_trials,3); % [ADD, MTFA, Tau]
results_CUSUM_D = zeros(n_trials,3); % [ADD, MTFA, Tau]
results_FILTER_B = zeros(n_trials,3); % [ADD, PFA, Tau]
results_FILTER_D = zeros(n_trials,3); % [ADD, PFA, Tau]

% Make a progress bar to see the simulation progress
u = waitbar(0, 'Simulation Progress');

for i = [1:n_trials]
   ADD_CUSUM_B = 0; MTFA_CUSUM_B = 0; tau_CUSUM_B = 0;
   ADD_CUSUM_D = 0; MTFA_CUSUM_D = 0; tau_CUSUM_D = 0;
   ADD_FILTER_B = 0; PFA_FILTER_B = 0; tau_FILTER_B = 0;
   ADD_FILTER_D = 0; PFA_FILTER_D = 0; tau_FILTER_D = 0;
   for j = [1:n_scenarios]
        %% RUN A BAYESIAN SCENARIO FIRST WITH A RANDOMLY GENERATED CHANGEPOINT
        %% Generate the Bayesian scenario

        % State sequence, observations and changepoint
        [X_B, y_B, nu_B] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
            means(i,:),var_affected);

        %% Fetch the algorithm's test statistics

        % ==== BAYESIAN ====
        % CUSUM
        % Set threshold:
        h_cusum = 8;
        [ADD,MTFA,tau,~] = CUSUMScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_B, nu_B, h_cusum);
         
        % Append results to SNR trial increment
        ADD_CUSUM_B = ADD_CUSUM_B + ADD;
        MTFA_CUSUM_B = MTFA_CUSUM_B + MTFA;
        tau_CUSUM_B = tau_CUSUM_B + tau;
        
        % HMM Filter
        % Set threshold:
        h_bay = 0.99;
        [ADD,PFA,tau,~] = BayesianScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_B, nu_B, h_bay);
        
        % Append results to SNR trial increment
        ADD_FILTER_B = ADD_FILTER_B + ADD;
        PFA_FILTER_B = PFA_FILTER_B + PFA;
        tau_FILTER_B = tau_FILTER_B + tau;

        %% GENERATED A DETERMINISTIC CHANGEPOINT SCENARIO

        % State sequence, observations and changepoint
        [X_D, y_D, nu_D] = simulate_deterministic_scenario(mean_unaffected,var_unaffected, ...
            means(i,:),var_affected);

        %% Fetch the deterministic algorithm's performance
        % ==== DETERMINISTIC ====
        % Assumes a rho of 1e-4
        % CUSUM
        [ADD,MTFA,tau,~] = CUSUMScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_D, nu_D, h_cusum);
        
        % Append results to SNR trial increment
        ADD_CUSUM_D = ADD_CUSUM_D + ADD;
        MTFA_CUSUM_D = MTFA_CUSUM_D + MTFA;
        tau_CUSUM_D = tau_CUSUM_D + tau;
        
        % HMM Filter
        [ADD,PFA,tau,~] = BayesianScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_D, nu_D, h_bay);
        
        % Append results to SNR trial increment
        ADD_FILTER_D = ADD_FILTER_D + ADD;
        PFA_FILTER_D = PFA_FILTER_D + PFA;
        tau_FILTER_D = tau_FILTER_D + tau;
        
        % Update the progress bar
        waitbar(((i*n_scenarios) + j)/(n_trials*n_scenarios));
   end
   
   % Take the means and insert into the main results matrix
   results_CUSUM_B(i,:) = [ADD_CUSUM_B MTFA_CUSUM_B tau_CUSUM_B] ./ n_scenarios;
   results_CUSUM_D(i,:) = [ADD_CUSUM_D MTFA_CUSUM_D tau_CUSUM_D] ./ n_scenarios;
   results_FILTER_B(i,:) = [ADD_FILTER_B PFA_FILTER_B tau_FILTER_B] ./ n_scenarios;
   results_FILTER_D(i,:) = [ADD_FILTER_D PFA_FILTER_D tau_FILTER_D] ./ n_scenarios;
end

% cleanup
clearvars ADD MTFA tau S M_hat PFA

%% Display results

% Stopping Time
figure

% Random changepoint
subplot(1,2,1)
hold on
colororder({'b','r'})
yyaxis left
plot(SNR, results_FILTER_B(:,1),'b*') % ADD
ylabel('HMM Filter Delay','Interpreter','Latex')
yyaxis right
plot(SNR, results_CUSUM_B(:,1),'r.') % PFA
ylabel('CUSUM Delay','Interpreter','Latex')
title('\textbf{Random}', ...
    'Interpreter', 'Latex')
xlabel('SNR (dB)')
set(gca, 'color', [0 0.07 0.1 0.2])

subplot(1,2,2)
hold on
colororder({'b','r'})
yyaxis left
plot(SNR, results_FILTER_D(:,1),'b*') % ADD
ylabel('HMM Filter Delay','Interpreter','Latex')
yyaxis right
plot(SNR, results_CUSUM_D(:,1),'r.') % PFA
ylabel('CUSUM Delay','Interpreter','Latex')
title('\textbf{Deterministic}', ...
    'Interpreter', 'Latex')
xlabel('SNR (dB)')
set(gca, 'color', [0 0.07 0.1 0.2])