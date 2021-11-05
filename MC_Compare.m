% This script is used to test the performance of the Bayesian and CUSUM
% techniques with a specific set of means and variances iaw a SNR bound
close all
clear
clc

% SET ENVIRONMENT SETTINGS
do_random_mean = 0;
random_mean_modifier = [ 1 1 1 ];

% Add the paths requried to run
addpath('DatasetGen')
addpath('Simulation')

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
SNR_min = -15;
SNR_max = -10;
n_trials = 20 + 1;

% Form the SNR vector
SNR = linspace(SNR_min, SNR_max, n_trials);

means = zeros(n_trials,length(mean_unaffected));
% Determine the new means to test with depending on the desired SNR
for i = [1:n_trials]
    means(i,:) = var_unaffected .* db2pow(SNR(i)) + mean_unaffected;
end

%% Begin iterating around each SNR and running a number of trials
n_scenarios = 200;

% Store the results in matrices
results_CUSUM_R = zeros(n_trials,3); % [ADD, MTFA, Tau]
results_CUSUM_D = zeros(n_trials,3); % [ADD, MTFA, Tau]
results_FILTER_R = zeros(n_trials,3); % [ADD, PFA, Tau]
results_FILTER_D = zeros(n_trials,3); % [ADD, PFA, Tau]
results_MC = zeros(n_trials,4); % [CUSUM_B, FILTER_B, CUSUM_D, FILTER_D]
results_FA = zeros(n_trials,4); % [CUSUM_B, FILTER_B, CUSUM_D, FILTER_D]

% Make a progress bar to see the simulation progress
u = waitbar(0, 'Simulation Progress');

% Thresholds
%h_B = [0.999769140523478,0.999969578001201,0.999795850507762,0.999816372444712,0.999644934675811,0.999867795754736,0.999001136230981,0.999263866186132,0.999897176211415,0.999166297478432,0.998416738053169];
%h_C = [13.0044624250654,13.5233154436896,13.2375955343193,13.2065419692670,13.6886603565447,15.1927294870466,13.6174437354134,13.9509832102166,16.9185020257958,17.5272155064804,15.5988191752556];


for i = [1:n_trials]
   ADD_CUSUM_R = 0; MTFA_CUSUM_R = 0; tau_CUSUM_R = 0;
   ADD_CUSUM_D = 0; MTFA_CUSUM_D = 0; tau_CUSUM_D = 0;
   ADD_FILTER_R = 0; PFA_FILTER_R = 0; tau_FILTER_R = 0;
   ADD_FILTER_D = 0; PFA_FILTER_D = 0; tau_FILTER_D = 0;
   MC_R_CUSUM = 0; MC_D_CUSUM = 0; MC_R_FILTER = 0; MC_D_FILTER = 0;
   FA_R_CUSUM = 0; FA_D_CUSUM = 0; FA_R_FILTER = 0; FA_D_FILTER = 0;
   for j = [1:n_scenarios]
        %% RUN A BAYESIAN SCENARIO FIRST WITH A RANDOMLY GENERATED CHANGEPOINT
        %% Generate the Bayesian scenario

        % State sequence, observations and changepoint
        [X_B, y_B, nu_B] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
            means(i,:),var_affected);

        %% Fetch the algorithm's test statistics

        % ==== RANDOM ====
        % CUSUM
        % Set threshold:
        h_cusum = 13.178396938960308;
        [ADD,MTFA,tau,~] = CUSUMScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_B, nu_B, h_cusum);
         
        % Append results to SNR trial increment
        if ADD == -1 || MTFA == -1 || tau == -1
            MC_R_CUSUM = MC_R_CUSUM + 1;
        elseif tau < nu_B
            FA_R_CUSUM = FA_R_CUSUM + 1;
        else
            ADD_CUSUM_R = ADD_CUSUM_R + ADD;
            MTFA_CUSUM_R = MTFA_CUSUM_R + MTFA;
            tau_CUSUM_R = tau_CUSUM_R + tau;
        end
        
        % HMM Filter
        % Set threshold:
        h_bay = 0.999701274462108;
        [ADD,PFA,tau,~] = BayesianScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_B, nu_B, h_bay);
        
        % Append results to SNR trial increment
        if ADD == -1 || PFA == -1 || tau == -1
            MC_R_FILTER = MC_R_FILTER + 1;
        elseif tau < nu_B
            FA_R_FILTER = FA_R_FILTER + 1;
        else
            ADD_FILTER_R = ADD_FILTER_R + ADD;
            PFA_FILTER_R = PFA_FILTER_R + PFA;
            tau_FILTER_R = tau_FILTER_R + tau;
        end

        %% GENERATED A DETERMINISTIC CHANGEPOINT SCENARIO

        % State sequence, observations and changepoint
        [X_D, y_D, nu_D] = simulate_deterministic_scenario(mean_unaffected,var_unaffected, ...
            means(i,:),var_affected, 1878);

        %% Fetch the deterministic algorithm's performance
        % ==== DETERMINISTIC ====

        % CUSUM
        h_cusum = 12.806962979038087;
        [ADD,MTFA,tau,~] = CUSUMScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_D, nu_D, h_cusum);
        
        % Append results to SNR trial increment
        if ADD == -1 || MTFA == -1 || tau == -1
            MC_D_CUSUM = MC_D_CUSUM + 1;
        elseif tau < nu_D
            FA_D_CUSUM = FA_D_CUSUM + 1;
        else
            ADD_CUSUM_D = ADD_CUSUM_D + ADD;
            MTFA_CUSUM_D = MTFA_CUSUM_D + MTFA;
            tau_CUSUM_D = tau_CUSUM_D + tau;  
        end
        
        % HMM Filter
        h_bay = 0.999623011766654;
        [ADD,PFA,tau,~] = BayesianScenario(mean_unaffected, ...
            var_unaffected, means(i,:), var_affected, y_D, nu_D, h_bay);
        
        % Append results to SNR trial increment
        if ADD == -1 || PFA == -1 || tau == -1
            MC_D_FILTER = MC_D_FILTER + 1;
        elseif tau < nu_D
            FA_D_FILTER = FA_D_FILTER + 1;
        else
            ADD_FILTER_D = ADD_FILTER_D + ADD;
            PFA_FILTER_D = PFA_FILTER_D + PFA;
            tau_FILTER_D = tau_FILTER_D + tau; 
        end
        
        % Update the progress bar
        waitbar((((i-1)*n_scenarios) + j)/(n_trials*n_scenarios));
   end
   
   % Take the means and insert into the main results matrix
   results_CUSUM_R(i,:) = [ADD_CUSUM_R MTFA_CUSUM_R tau_CUSUM_R] ./ n_scenarios;
   results_CUSUM_D(i,:) = [ADD_CUSUM_D MTFA_CUSUM_D tau_CUSUM_D] ./ n_scenarios;
   results_FILTER_R(i,:) = [ADD_FILTER_R PFA_FILTER_R tau_FILTER_R] ./ n_scenarios;
   results_FILTER_D(i,:) = [ADD_FILTER_D PFA_FILTER_D tau_FILTER_D] ./ n_scenarios;
   results_MC(i,:) = [MC_R_CUSUM MC_R_FILTER MC_D_CUSUM MC_D_FILTER];
   results_FA(i,:) = [FA_R_CUSUM FA_R_FILTER FA_D_CUSUM FA_D_FILTER];
end

close(u)

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
plot(SNR, results_FILTER_R(:,1),'b*') % ADD
ylabel('HMM Filter Delay','Interpreter','Latex')
ylim([0 1.05 * max([results_FILTER_R(:,1) results_CUSUM_R(:,1)],[],'ALL')])
yyaxis right
plot(SNR, results_CUSUM_R(:,1),'r.','MarkerSize',12) % PFA
ylabel('CUSUM Delay','Interpreter','Latex')
ylim([0 1.05 * max([results_FILTER_R(:,1) results_CUSUM_R(:,1)],[],'ALL')])
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
ylim([0 1.05 * max([results_FILTER_D(:,1) results_CUSUM_D(:,1)],[],'ALL')])
yyaxis right
plot(SNR, results_CUSUM_D(:,1),'r.','MarkerSize',12) % PFA
ylabel('CUSUM Delay','Interpreter','Latex')
ylim([0 1.05 * max([results_FILTER_D(:,1) results_CUSUM_D(:,1)],[],'ALL')])
title('\textbf{Deterministic}', ...
    'Interpreter', 'Latex')
xlabel('SNR (dB)')
set(gca, 'color', [0 0.07 0.1 0.2])