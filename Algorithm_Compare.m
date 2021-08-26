% This script is used to test the performance of the Bayesian and CUSUM
% techniques with a specific set of means and variances
close all
clear
clc

% SET ENVIRONMENT SETTINGS
do_random_mean = 1;
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

var_affected = [1 1 1];

% Print simulation parameters
disp("=============================================")
disp("Running comparison simulation with parameters:")
disp(" ")
disp("Means (Unaffected) = [" + num2str(mean_unaffected) + "],")
disp("Variances (Unaffected) = [" + num2str(var_unaffected) + "],")
disp("Means (Affected) = [" + num2str(mean_affected) + "],")
disp("Variances (Affected) = [" + num2str(var_affected) + "].")
disp("=============================================")
disp(" ")

%% RUN A BAYESIAN SCENARIO FIRST WITH A RANDOMLY GENERATED CHANGEPOINT
%% Generate the Bayesian scenario

% State sequence, observations and changepoint
[X_B, y_B, nu_B] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected);

%% Fetch the algorithm's test statistics

% ==== BAYESIAN ====
% CUSUM
% Set threshold:
h_cusum = 5;
[ADD_CUSUM_B,PFA_CUSUM_B,tau_CUSUM_B,S_B] = CUSUMScenario(mean_unaffected, ...
    var_unaffected, mean_affected, var_affected, y_B, nu_B, h_cusum);

% HMM Filter
% Set threshold:
h_bay = 0.99;
[ADD_HMM_B,MTFA_HMM,tau_HMM_B,M_hat_B] = BayesianScenario(mean_unaffected, ...
    var_unaffected, mean_affected, var_affected, y_B, nu_B, h_bay);

%% GENERATED A DETERMINISTIC CHANGEPOINT SCENARIO

% State sequence, observations and changepoint
[X_D, y_D, nu_D] = simulate_deterministic_scenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected, nu_B);

% ==== DETERMINISTIC ====
% Deterministic CUSUM
h_cusum_D = 15;
[ADD_CUSUM_D,D_PFA_CUSUM_B,tau_CUSUM_D,S_D] = Deterministic_CUSUMScenario(mean_unaffected, ...
    var_unaffected, mean_affected, var_affected, y_D, h_cusum_D, nu_D);

% Deterministic HMM Filter
% [ADD_CUSUM_D,D_PFA_CUSUM_B,tau_CUSUM_D,S_D] = Deterministic_CUSUMScenario(mean_unaffected, ...
%     var_unaffected, mean_affected, var_affected, X_D, y_D, nu_D, h_cusum);
%% Plot the findings

% Number of samples
n_samples = length(X_B);

figure

% Bayesian Plots

% CUSUM Plot
subplot(2,2,1)
hold on
plot([1:n_samples],S_B) % Plot likelihood of system in post-change
plot(nu_B,S_B(nu_B),'go') % Plot the change-point
plot(tau_CUSUM_B,S_B(tau_CUSUM_B),'ro') % Plot the stopping time
yline(h_cusum,'m--') % Plot the detection threshold
hold off

set(gca, 'color', [0 0.07 0.1 0.2])
set(gca, 'YScale', 'log')
title('Bayesian CUSUM Algorithm $$S_k$$','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$S_k$$','Interpreter','Latex')
xlim([0 n_samples])
ylim([0 5000])

% HMM Filter Plot
subplot(2,2,2)
hold on
plot([1:n_samples],M_hat_B(2,:)) % Plot likelihood of system in post-change
plot(nu_B,M_hat_B(2,nu_B),'go') % Plot the change-point
plot(tau_HMM_B,M_hat_B(2,tau_HMM_B),'ro') % Plot the stopping time
yline(h_bay,'m--') % Plot the detection threshold
hold off

set(gca, 'color', [0 0.07 0.1 0.2])
title('Bayesian HMM Algorithm $$\hat{M}_k^2$$','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$M_k^2$$','Interpreter','Latex')
xlim([0 n_samples])
ylim([-0.1 1.1])


% Deterministic Plots

% CUSUM plot
subplot(2,2,3)
hold on
plot([1:n_samples],S_D) % Plot likelihood of system in post-change
plot(nu_D,S_D(nu_D),'go') % Plot the change-point
plot(tau_CUSUM_D,S_D(tau_CUSUM_D),'ro') % Plot the stopping time
yline(h_cusum,'m--') % Plot the detection threshold
hold off

set(gca, 'color', [0 0.07 0.1 0.2])
set(gca, 'YScale', 'log')
title('Determinstic CUSUM Algorithm $$S_k$$','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$S_k$$','Interpreter','Latex')
xlim([0 n_samples])
ylim([0 5000])

%% Print results

disp(" ")
disp("RESULTS:")
disp("It took the Bayesian CUSUM algorithm " + num2str(ADD_CUSUM_B) + ... 
    " samples to detect the change event");
disp("It took the Bayesian HMM algorithm " + num2str(ADD_HMM_B) + ...
    " samples to detect the change event");
disp("It took the Deterministic CUSUM algorithm " + num2str(ADD_CUSUM_D) + ...
    " samples to detect the change event");
