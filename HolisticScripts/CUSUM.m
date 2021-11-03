%% Environment Clearing
clc
close all
clearvars

addpath('..\Helper');
addpath('..\DatasetGen');
addpath('..\Simulation');

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

% Cleanup
clearvars trans_beta

%% Generate the state randoms and changepoint

% Check whether the changepoint is already supplied
if ~exist('nu', 'var')
    % Generate a random changepoint between 0 and n_samples / 2
    nu = randi([1 n_samples]);
end

% Check if nu is invalid
while nu >= n_samples
    disp("Regenerating changepoint of deterministic scenario")
    nu = randi([1 n_samples]);
end

% Define changepoint as the number of samples
%nu = n_samples;

% Initialise the state sequence with ones until the changepoint to indicate
% the system is in the pre-change state
X = ones(1, n_samples);

if nu < n_samples
    % Simulate the rest of the markov chain in the post-change state space
    % using the post-change transition matrices
    X(nu:end) = simulate(A_beta, n_samples - nu) + 1;
end

%% Define the probabilities of specific events

% There is no pre-change state dependence, such that rho is constant for
% all states in the state spaces

% Rho is the probability that an object will enter the network
rho = 5e-4;

% The probability that the network will tranisition from state alpha to
% beta
A_nu = rho * pi_k;

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

% Fetch the transition points from the state sequence X
trans = getTransitionIndices(X);

%% Define the Distribution Parameters

% Assume the sensors observation measurements are i.i.d.

% Each sensor's default distribution without being affected by a target
% will be a standard normal distribution with sigma = 1 and mean = 0
mean_unaffected = [1,2,3]; %= [1 2 3];
var_unaffected = [1 1 1];

% Each distribution parameter can be accessed with the node's index
%variances = dist_dev + 0.5*rand(1,n_states);
mean_affected = [2,3,4];%= [2 3 4];
var_affected = [1 1 1];

%% Generate randomly distributed values for each sensing node

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors,n_samples);

% Populate the y matrix with samples taken from a Gaussian distribution
% with mean and variances per state as defined in means and vars
for i = [1:n_samples]
    % Check whether we are in pre or post-change
    if(i < nu) || nu == n_samples
        % Generate an unaffected distribution sample
        y(:,i) = sqrt(var_unaffected).' .* randn(n_sensors,1) + ...
            mean_unaffected.';
    else
        % Define the var and mean to generate data with at this sample time
        means = mean_unaffected; means(X(i)-1) = mean_affected(X(i)-1);
        vars = var_unaffected; vars(X(i)-1) = var_affected(X(i)-1);

        % Generate the data using a scaled randn value
        y(:,i) = sqrt(vars).' .* randn(n_sensors,1) + means.';
    end
end

% Plot the observation data
plotObservationData(n_sensors,trans,y,nu,mean_unaffected);

%% Hidden Markov Model Filter

% To estimate the current state of the HMM from the observation densities,
% the HMM filter is calculated in O(n^2)
 
Z_hat = zeros(n_sensors,n_samples);
Z_hat(:,1) = pi_k.';
 
% Set the test statistic to start at node 1, which is the pre-change state
% Z_hat(:,1) = [1 zeros(1,n_sensors)].';
 
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

% Define a probability threshold to test for
h = 8;

% Form a set of k values
k_h = [1:n_samples];

% Index the set such that it contains entries of S > h
k_h = k_h(S > h);

% Determine the lower bound of the M_h set to determine the first location
% of when the test statistic exceeds the probability threshold
tau = min(k_h);

% Cleanup
clearvars k_h

%% Plot the stopping results

plotCUSUMResults(S, h, nu, tau);

%% Calculate performance parameters

% % Average Detection Delay
ADD = max(0,tau - nu);

% % Probability of False Alarm

% Take the assumption that the stopping time is the first time the CUSUM
% statistic has crossed the threshold since the algorithms execution
B_comp = -(mean(S(1:tau)) - tau) / (tau + 1);
PFA = (1 - B_comp) / (mean(S(1:tau)) - B_comp);

%% Cleanup
clearvars i j