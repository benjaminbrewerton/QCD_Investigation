function [ADD,MTFA,tau,S] = CUSUMScenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected, y, nu, h)
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

%% Define Transition Matrices for states

% Form the Discrete Time Markov Chains for each sensor space

% Since in state alpha, a singular non-sensed state is modelled:
A_alpha = dtmc([1]);

% Define a transition modifier which models what the probability is of a
% detected object moving outside of the current state
P_move = 0.01;

% In state beta, there are n_sensors nodes, so to create a random
% transition matrix
% Uncomment the line below for a fully randomised state transition matrix
%A_beta = mcmix(n_sensors);

% This beta DTMC is formed from the P_move variable on the diagonal and all
% other values are P_move/2
trans_beta = diag(repelem(1 - P_move, n_sensors));
trans_beta(trans_beta == 0) = P_move / (n_sensors - 1);
A_beta = dtmc(trans_beta);

% For transition state nu, assume that the probability for the state to
% tranisition to the states in space beta is equal
pi_k = ones(1,n_sensors) ./ n_sensors;

%% Hidden Markov Model Filter

% To estimate the current state of the HMM from the observation densities,
% the HMM filter is calculated in O(n^2)
 
Z_hat = zeros(n_sensors,n_samples);
Z_hat(:,1) = pi_k.';
 
% Initialise an array S, which contains the CUSUM test statistic
S = zeros(1,n_samples);

for i = [2:n_samples]
    % Get the observation vector for a time sample k
    cur_obs = y(:,i);
    
    % Define a square matrix B whose values whose rows represent the
    % probability of a singular observation being from the set of all
    % densities in the system
    B = zeros(n_sensors);
    
    for j = [1:n_sensors]
        % Initialise the means and variances for each element
        cur_vars = var_unaffected;
        cur_means = mean_unaffected;
 
        % Modify the mean and dists in position j to reflect the mean 
        % of the affected distributions
        cur_means(j) = mean_affected(j);
        cur_vars(j) = var_affected(j);
 
        % Populate with the affected distribution
        B(:,j) = (1./sqrt(2*pi*cur_vars)).' .* ...
            exp(-(cur_obs - cur_means.').^2 ./ (2*cur_vars.'));
    end
 
    % Calculate the B matrix which is the diagonal of the PDF values at
    % each observation value for a sample k (i)
    B = diag(prod(B,1));
    
    % Determine the pre-change density likelihood
    P_alpha = prod((1./sqrt(2*pi*var_unaffected)).' .* ... 
        exp(-(cur_obs - mean_unaffected.').^2 ./ (2*var_unaffected.')));
    
    % Get the previously calculated test statistic
    Z_prev = Z_hat(:,i-1);
    
    % Define the new Z test statistic
    Z_new = B * A_beta.P.' * Z_prev; % Big A matrix
    
    % Calculate the normalistation factor
    N = 1 / sum(Z_new);
    
    % Insert the new state estimate
    Z_hat(:,i) = N * Z_new;
    
    % Evaluate the test statistic using the CUSUM algorithm under Lorden's
    % criteria
    S(i) = max(0, S(i-1) + log(1/N) - log(P_alpha));
end

%% Infimum Bound Stopping Time

% Check if threshold exists
if ~exist('h','var')
    h = 8;
end

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

% Mean time to false alarm
%PFA = 1 - M_hat(2,tau);
MTFA = n_samples;

if isempty(tau)
    tau = -1;
end
if isempty(ADD)
    ADD = -1;
end
end