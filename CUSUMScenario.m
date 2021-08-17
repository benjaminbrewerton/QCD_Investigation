function [ADD,PFA,tau,S] = CUSUMScenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected, X, y, nu, h)
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

%% Define the probabilities of specific events

% There is no pre-change state dependence, such that rho is constant for
% all states in the state spaces

% Rho is the probability that an object will enter the network
rho = 5e-4;

% The probability that the network will tranisition from state alpha to
% beta
A_nu = rho * pi_k;

%% Generate the state randoms and changepoint

% Since before the changepoint, the state variables in space alpha are not
% relevant to the problem and do not have to be simulated.

% Calculate A to determine the state variables X
% Use the definition where rho is constant for all states
A = dtmc([(1-rho)*A_alpha.P A_nu ; zeros(n_sensors,1) A_beta.P], ...
    'StateNames',["Alpha 1" "Beta 1" "Beta 2" "Beta 3"]);

%% Determine the transition points

% Define a new matrix, e, which will hold a 1 in the row where the current
% position in the state vector is
e = zeros(n_states,n_samples);

% Loop around the state vector and place a 1 in the e vector in the
% appropriate column
for i = [1:n_samples]
   e(X(i),i) = 1; 
end

%% Hidden Markov Model Filter

% To estimate the current state of the HMM from the observation densities,
% the HMM filter is calculated in O(n^2)
 
Z_hat = zeros(n_states,n_samples);
 
% Set the test statistic to start at node 1, which is the pre-change state
Z_hat(:,1) = [1 zeros(1,n_sensors)].';
 
% Initialise an array S, which contains the CUSUM test statistic
S = zeros(1,n_samples);
 
% Transpose the A matrix to align with literature definition
AT = A.P.';
 
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
    B = diag(prod(B,1));
    
    % Get the previously calculated test statistic
    Z_prev = Z_hat(:,i-1);
    
    % Define the new Z test statistic
    Z_new = B * AT * Z_prev; % Big A matrix
    
    % Calculate the normalistation factor
    N = 1 / sum(Z_new);
    
    % Insert the new state estimate
    Z_hat(:,i) = N * Z_new;
    
    % Evaluate the test statistic using the CUSUM algorithm under Lorden's
    % criteria
    S(i) = max(0, S(i-1) + log(1/N) - log(B(1,1)));
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