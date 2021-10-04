%% Environment Clearing
clc
close all
clearvars

addpath('..');

%% Begin definition of network variables

% Number of pre-change states
n_states_mu = 1;

% Number of states in the post-change state
n_states_nu = 3;

% Number of sensors
n_sensors = max(n_states_mu, n_states_nu);

% Total number of states
n_states = n_states_mu + n_states_nu;

% Number of random samples
n_samples = 1e4;

%% Define Transition Matrices for states

% Form the Discrete Time Markov Chains for each change space

% Define a transition modifier which models what the probability is of a
% detected object moving outside of the current state
% Assume that the probability of transitioning to an adjacent state is
% equivalent for both spaces.
P_move = 0.01;

% Form the DTMCs for each change space using the definition from ISD papers
trans = diag(repelem(1 - P_move, n_states_mu));
trans(trans == 0) = P_move / (n_states_mu - 1);
A_alpha = dtmc(trans);
% Assume both the change space DTMCs operate with the same probabilities of
% transition
trans = diag(repelem(1 - P_move, n_states_nu));
trans(trans == 0) = P_move / (n_states_nu - 1);
A_beta = dtmc(trans);

% Also assume the probability for beginning in each state in the DTMCs is
% equal
pi_mu = ones(1,n_states_mu) ./ n_states_mu;
pi_nu = ones(1,n_states_nu) ./ n_states_nu;

% Cleanup
clearvars trans

%% Generate the state randoms and changepoint

% Generate a specific number of changepoints
n_changepoints = 10;

% Generate the checker pattern
checkers = ones(1,n_changepoints);
for i = [1:n_changepoints]
    if mod(i,2) == 0
        checkers(i) = -1;
    end
end

% Check whether the changepoint is already supplied
if ~exist('lambda', 'var')
    % Generate a random changepoint between 0 and n_samples / 2
    lambda = sort(randi([1 n_samples],1,n_changepoints)) .* checkers;
end

% Ensure that the maximum changepoint value is n_samples-1
while max(lambda) >= n_samples
    lambda = sort(randi([1 n_samples],1,n_changepoints)) .* checkers;
end

% Define changepoint as the number of samples
%nu = [];

% Initialise the state sequence with ones until the changepoint to indicate
% the system is in the pre-change state
X = ones(1, n_samples);

% Loop around the changepoints and assign simulated values
for i = [1:length(lambda)]
    % Get the current changepoint
    lambda_cur = lambda(i);
    % Get the next changepoint
    if i ~= length(lambda)
        lambda_next = abs(lambda(i+1));
    else
        lambda_next = n_samples;
    end
    
    % Assign appropriate state values for when the changepoint is in the
    % post-change space
    if lambda_cur > 0
        X(lambda_cur:lambda_next) = simulate(A_beta, lambda_next-lambda_cur) + 1;
    end
end

% Print the generated changepoints
disp("Generated scenario with changepoints at: ")
for i = [1:length(lambda)]
   if i ~= length(lambda)
        disp([char(9) num2str(abs(lambda(i))) ',']);
   else
        disp([char(9) num2str(abs(lambda(i))) '.']);
   end
end
%% Determine the transition points for plotting

% Fetch the transition points from the state sequence X
trans = getTransitionIndices(X);

%% Define the Distribution Parameters

% Assume the sensors observation measurements are i.i.d.

% Each sensor's default distribution without being affected by a target
% will be a standard normal distribution with sigma = 1 and mean = 0
mean_mu = [1,2,3]; %= [1 2 3];
var_mu = [1 1 1];

% Each distribution parameter can be accessed with the node's index
%variances = dist_dev + 0.5*rand(1,n_states);
mean_nu = [2,3,4];%= [2 3 4];
var_nu = [1 1 1];

%% Generate randomly distributed values for each sensing node

% Here we assume the number of sensors is equal in the pre-change and
% post-change scenarios

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors, n_samples);

% Populate the y matrix with samples taken from a Gaussian distribution
% with mean and variances per state as defined in means and vars
for i = [1:n_samples]
    % Check whether we are in pre-change or post-change states
    if(X(i) <= n_states_mu)
        % Generate an unaffected distribution sample
        y(:,i) = sqrt(var_mu).' .* randn(n_sensors,1) + ...
            mean_mu.';
    else
        % Define the var and mean to generate data with at this sample time
        means = mean_mu; means(X(i)-1) = mean_nu(X(i)-1);
        vars = var_mu; vars(X(i)-1) = var_nu(X(i)-1);

        % Generate the data using a scaled randn value
        y(:,i) = sqrt(vars).' .* randn(n_sensors,1) + means.';
    end
end

%% Plot the observation data

plotObservationData(n_sensors, trans, y, lambda, mean_mu)

%% Hidden Markov Model Filter

% To estimate the current state of the HMM from the observation densities,
% the HMM filter is calculated in O(n^2)

% Initialise the test statistic to be equal probability of occurrence
Z_hat = zeros(n_sensors,n_samples);
Z_hat(:,1) = pi_nu.';
 
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
        cur_vars = var_mu;
        cur_means = mean_mu;
 
        % Modify the mean and dists in position j to reflect the mean 
        % of the affected distributions
        cur_means(j) = mean_nu(j);
        cur_vars(j) = var_nu(j);
 
        % Populate with the affected distribution
        B(:,j) = (1./sqrt(2*pi*cur_vars)).' .* ...
            exp(-(cur_obs - cur_means.').^2 ./ (2*cur_vars.'));
    end
 
    % Calculate the B matrix which is the diagonal of the PDF values at
    % each observation value for a sample k (i)
    B = diag(prod(B,1));
    
    % Determine the pre-change density likelihood
    P_alpha = prod((1./sqrt(2*pi*var_mu)).' .* ... 
        exp(-(cur_obs - mean_mu.').^2 ./ (2*var_mu.')));
    
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

% Cleanup
clearvars Z_prev Z_new Z_ins B_cur

%% Define the Mode Process Vector

% Initialise an empty 2 x n_samples matrix to store the mode process
% variable
M_hat = zeros(2, n_samples);

% Use indexing to calculate M_k = M(Z_k)
M_hat(1,:) = Z_hat(1,:);
M_hat(2,:) = 1 - Z_hat(1,:);

%% Infimum Bound Stopping Time

% Define the threshold for the CUSUM statistic
h = 8;

% Initialise a vector which stores the threshold crossing indices. Store as
% -1 if no changepoint is detected within the common boundary.
tau = ones(1,length(lambda)) .* -1;

% Loop around the changepoints and determine the index of the first sample
% which exceeds the specified thresholds
for i = [1:length(lambda)]
   % Get the current changepoint start index
   lambda_cur = lambda(i);
   
   % Get the next changepoint index
   if i ~= length(lambda)
       lambda_next = lambda(i + 1);
   else
       lambda_next = n_samples;
   end
  
   % Begin iterating around each sample, looking for the point of threshold
   % crossing
   j = abs(lambda_cur);
   if lambda_cur > 0 % Pre-change to post-change case
    while j < abs(lambda_next) && ...
            S(j) < S(abs(lambda_cur)) + h
        j = j + 1; % Increment search index
    end
   else % Post-change to pre-change
    while j < abs(lambda_next) && ...
            S(j) > S(abs(lambda_cur)) - h
        j = j + 1; % Increment search index
    end
   end
    
    % Input the changepoint determined into the primary set
    if j ~= abs(lambda_next)
        tau(i) = j;
    end
end

% Cleanup
clearvars M_stat

%% Plot the stopping results

plotCUSUMResults(S, [], lambda, tau);

%% Calculate performance parameters

% Average Detection Delay
DD = abs(abs(lambda) - tau);
DD(tau == -1) = -1;

ADD = mean(abs(abs(lambda) - tau));

% Probability of False Alarm
%PFA = 1 - M_hat(2,tau);

%% Cleanup
clearvars i j