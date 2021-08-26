function [ADD,PFA,tau,S] = Deterministic_CUSUMScenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected, y, h, nu)
% Produces the Average Delay Detection and Probability of False Alarm from
% a set of inputted means and variances (unaffected and affected
% respectively in order). This function assumes a setup of 3 sensors in the
% network.

%% Begin definition of network variables

% Number of pre-change states
n_states_pre_change = 1;

% Number of sensors for post-change modelling
n_sensors = 3;

% Total number of states
n_states = n_states_pre_change + n_sensors;

% Number of random samples
n_samples = 1e4;

%% Hidden Markov Model Filter
 
% Initialise an array S, which contains the CUSUM test statistic
S = zeros(1,n_samples);

for i = [2:n_samples]
    % Get the observation vector for a time sample k
    cur_obs = y(:,i);
    
    % Define a square matrix B whose values whose rows represent the
    % probability of a singular observation being from the set of all
    % densities in the system
    B = zeros(n_sensors,n_states);
    
    for j = [1:n_states]
        % Initialise the means and variances for each element
        cur_vars = var_unaffected;
        cur_means = mean_unaffected;
 
        if j ~= 1
        % Modify the mean and dists in position j to reflect the mean 
        % of the affected distributions
        cur_means(j-1) = mean_affected(j-1);
        cur_vars(j-1) = var_affected(j-1);
        end
 
        % Populate with the affected distribution
        B(:,j) = exp(-(cur_obs - cur_means.').^2 ./ (2*cur_vars.'));
    end
 
    % Calculate the B matrix which is the diagonal of the PDF values at
    % each observation value for a sample k (i)
    B = prod(B,1);
    
    % Evaluate the test statistic using the CUSUM algorithm under Lorden's
    % criteria
    S(i) = max(0, S(i-1) + log(mean(B(2:end))) ...
        - log(B(1,1)));
end

%% Infimum Bound Stopping Time

% Form a set of k values
k_h = [1:n_samples];

% Index the set such that it contains entries of S > h
k_h = k_h(S > h);

% Determine the lower bound of the M_h set to determine the first location
% of when the test statistic exceeds the probability threshold
tau = min(k_h);

% Cleanup
clearvars k_h

%% Calculate performance parameters

% Average Detection Delay
ADD = max(0,tau - nu);

% Probability of False Alarm
%PFA = 1 - M_hat(2,tau);
PFA = 0;

end