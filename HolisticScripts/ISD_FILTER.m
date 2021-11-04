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

%% Bayesian Probability Space Assumptions

% These variables are used for Bayesian changepoint modelling.

% Rho is the probability that an object will transition between state
% spaces.
rho = 2.5e-3;

% The probability that the network will tranisition from state alpha to
% beta
A_mu = rho .* pi_mu;
A_nu = rho .* pi_nu;

%% Generate the state randoms and changepoint

% Define the A from literature
A_lit = [(1-rho).*A_alpha.P (1/n_states_nu).*repelem(A_mu,n_states_mu,n_states_nu) ; ...
    repelem(A_nu.',n_states_mu,1) (1-rho .* (1/n_states_nu)).*A_beta.P];

% Determine the large A matrix which governs 
A = dtmc(A_lit.','StateNames',["Mu 1","Nu 1", "Nu 2", "Nu 3"]);

% Display the transition probabilities of the overall system
figure
graphplot(A,'ColorEdges',true)
title('Holistic ISD System Transition Probabilities')

% Simulate the augmented HMM established in matrix A.
% Start from node 1 in the pre-change space to simulate a sensor network
% that has not yet detected a target.
X = simulate(A, n_samples - 1,'X0',[1 zeros(1,n_states-1)]);

% Determine the points at which the HMM transitions between state spaces
lambda = zeros(1,n_samples);
lambda(X > n_states_mu) = 1; % Post-change
lambda(X <= n_states_mu) = 0; % Pre-change

% Get the differences in the changepoint vector to determine where these
% peaks occur. Append a 0 to the front of the vector to indicate the
% changepoint is detected on the sample of occurrence instead of before.
% Store as 1 (moving from mu to nu) and -1 (moving from nu to mu).
lambda = [0 diff(lambda)];

% Get the state indices of changepoints
lambda = [1:n_samples] .* lambda;
% Remove all zero indices
lambda(lambda == 0) = [];


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
n_sensors = max(n_states_mu,n_states_nu);

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

% Initialise the filter recursion by setting the initial conditions of the
% estimated state sequence to be the distribution of initial DTMC node. As 
% with previous formatting, each row will represent the test statistics for
% a particular sensor node whilst the columns represent the test statistics
% for each sensor node at a particular time instance, k
Z_hat = zeros(n_states,n_samples);

% Initialise the Z test to have the equal distribution pi as the first
% column entry
Z_hat(:,1) = [1 zeros(1,n_sensors)].';

for i = [2:n_samples]
    % Get the observation vector for a time sample k
    cur_obs = y(:,i);

    % Define a square matrix B_cur of B values whose rows represent the
    % probability of a singular observation being from the set of all
    % densities in the system
    B = zeros(n_sensors,n_states);
    
    for j = [1:n_states]
        % Initialise the means and variances for each element
        cur_vars = var_mu;
        cur_means = mean_mu;

        if j ~= 1
            % Modify the mean and dists in position j to reflect the mean 
            % of the affected distributions
            cur_means(j-1) = mean_nu(j-1);
            cur_vars(j-1) = var_nu(j-1);
        end

        % Populate with the affected distribution
        B(:,j) = (1./sqrt(2*pi*cur_vars)).' .* ...
            exp(-(cur_obs - cur_means.').^2 ./ (2*cur_vars.'));
    end
    
    % Calculate the B matrix which is the diagonal of the PDF values at
    % each observation value for a sample k (i)
    B = diag(prod(B,1));

    % Get the previous estimate
    Z_prev = Z_hat(:,i-1);

    % Calculate the new estimate
    Z_new = B * A_lit * Z_prev; % Big A matrix

    % Normalise the new estimate
    Z_ins = inv(sum(Z_new)) * Z_new;

    % Input the new estimate into the main Z test matrix
    Z_hat(:,i) = Z_ins;
end

% Cleanup
clearvars Z_prev Z_new Z_ins B_cur

%% Plot the test statistic results

plotTestStatistics(Z_hat, trans);

%% Alternate Z_k plot

plotTestAccuracy(Z_hat, trans, lambda);

%% Define the Mode Process Vector

% Initialise an empty 2 x n_samples matrix to store the mode process
% variable
M_hat = zeros(2, n_samples);

% Use indexing to calculate M_k = M(Z_k)
M_hat(1,:) = Z_hat(1,:);
M_hat(2,:) = 1 - Z_hat(1,:);

%% Infimum Bound Stopping Time

% Define a probability threshold to test for
h = 0.99;

% Form a set of k values
k_h = [1:n_samples];

% Get the mode statistic for when the system is in post-change
M_mu = M_hat(1,:);
M_nu = M_hat(2,:);

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
   if lambda_cur > 0 % Pre-change to post-change case
       M_stat = M_nu;
   else % Post-change to pre-change
       M_stat = M_mu;
   end
   
    j = abs(lambda_cur);
    while j < abs(lambda_next) && M_stat(j) <= h
        j = j + 1; % Increment search index
    end
    
    % Input the changepoint determined into the primary set
    if j ~= abs(lambda_next)
        tau(i) = j;
    end
end

% Cleanup
clearvars M_stat

%% Plot the stopping results

plotStoppingResults(lambda, tau, M_hat, h);

%% Calculate performance parameters

% Average Detection Delay
ADD = mean(abs(abs(lambda) - tau));

% Probability of False Alarm
%PFA = 1 - M_hat(2,tau);

%% Cleanup
clearvars i j