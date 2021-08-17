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
disp("Running simulation with parameters:")
disp(" ")
disp("Means (Unaffected) = [" + num2str(mean_unaffected) + "],")
disp("Variances (Unaffected) = [" + num2str(var_unaffected) + "],")
disp("Means (Affected) = [" + num2str(mean_affected) + "],")
disp("Variances (Affected) = [" + num2str(var_affected) + "].")
disp(" ")

%% Generate the scenario

% State sequence, observations and changepoint
[X, y, nu] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected);

%% Fetch the algorithm's test statistics

% CUSUM
% Set threshold:
h_cusum = 5;
[ADD_CUSUM,PFA_CUSUM,tau_CUSUM,S] = CUSUMScenario(mean_unaffected, ...
    var_unaffected, mean_affected, var_affected, X, y, nu, h_cusum);

% Bayesian
% Set threshold:
h_bay = 0.99;
[ADD_BAY,PFA_BAY,tau_BAY,M_hat] = BayesianScenario(mean_unaffected, ...
    var_unaffected, mean_affected, var_affected, X, y, nu, h_bay);

%% Plot the findings

% Number of samples
n_samples = length(X);

figure

% CUSUM Plot
subplot(1,2,1)
hold on
plot([1:n_samples],S) % Plot likelihood of system in post-change
plot(nu,S(nu),'go') % Plot the change-point
plot(tau_CUSUM,S(tau_CUSUM),'ro') % Plot the stopping time
yline(h_cusum,'m--') % Plot the detection threshold
hold off

set(gca, 'color', [0 0.07 0.1 0.2])
set(gca, 'YScale', 'log')
title('CUSUM Algorithm $$S_k$$','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$S_k$$','Interpreter','Latex')
% leg = legend('$$Z_k$$ -- $$\log\left(\frac{p_{\beta}\left(N^{-1}_{k|\lambda}\right)}{p_{\alpha}\left(N^{-1}_{k|\lambda}\right)}\right)$$',...
%     '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time', ...
%     '$$h$$ -- Threshold');
% leg.Interpreter = 'Latex';
% leg.Color = 'w';
xlim([0 n_samples])
ylim([0 5000])

% Bayesian Plot
subplot(1,2,2)
hold on
plot([1:n_samples],M_hat(2,:)) % Plot likelihood of system in post-change
plot(nu,M_hat(2,nu),'go') % Plot the change-point
plot(tau_BAY,M_hat(2,tau_BAY),'ro') % Plot the stopping time
yline(h_bay,'m--') % Plot the detection threshold
hold off

set(gca, 'color', [0 0.07 0.1 0.2])
title('Bayesian Algorithm $$\hat{M}_k^2$$','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$M_k^2$$','Interpreter','Latex')
% leg = legend('$$\hat{M}_k^2$$ -- $$P(X\in S_{\beta})$$',...
%     '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time', ...
%     '$$h$$ -- Threshold');
% leg.Interpreter = 'Latex';
% leg.Color = 'w';
xlim([0 n_samples])
ylim([-0.1 1.1])

%% Print results

disp(" ")
disp("RESULTS:")
disp("It took the CUSUM algorithm " + num2str(ADD_CUSUM) + ... 
    " samples to detect the change event");
disp("It took the Bayesian algorithm " + num2str(ADD_BAY) + ...
    " samples to detect the change event");