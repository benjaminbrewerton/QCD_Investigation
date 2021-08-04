%% Environment Clearing
clc
close all
clearvars

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

% Plot the Beta Sensor space transition probabilities
figure
graphplot(A_beta,'ColorEdges',true,'LabelEdges',true)
title('Post-Change Transition Probabilities')

% Cleanup
clearvars trans_beta

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

% Display the transition probabilities of the overall system
figure
graphplot(A,'ColorEdges',true)
title('Holistic System Transition Probabilities')

% Simulate the entire scenarios markov chain state transitions
% Assume that the when the space transitions from alpha to beta that the
% simulation begins at a random state with equal probability of initial
% state
% Initialise the sequence to begin at node 1
X = simulate(A, n_samples - 1,'X0',[1 zeros(1,n_sensors)]);

% Determine nu (changepoint) as the first time the sequence escapes from
% DTMC alpha and into DTMC beta
nu = length(X(X == 1)) + 1; % + 1 for being inclusive of the transition sample
% Print the result
disp(['It took ' num2str(nu) ...
    ' iterations to transition to the post-change state']);

% Check if the changepoint never occurs
if nu >= n_samples
   disp(["Error: Changepoint was never reached after " num2str(n_samples) ...
       " samples. Try again."]);
   quit
end
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
    if(i < nu)
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

% Initialise the filter recursion by setting the initial conditions of the
% estimated state sequence to be the distribution of initial DTMC node. As 
% with previous formatting, each row will represent the test statistics for
% a particular sensor node whilst the columns represent the test statistics
% for each sensor node at a particular time instance, k
Z_hat = zeros(n_states,n_samples);

% Initialise the Z test to have the equal distribution pi as the first
% column entry
Z_hat(:,1) = [1 zeros(1,n_sensors)].';

% Generate the transpose of the A matrix such that is inline with the
% defintion provided in literature
AT = A.P.';

for i = [2:n_samples]
    % Get the observation vector for a time sample k
    cur_obs = y(:,i);

    % Define a square matrix B_cur of B values whose rows represent the
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

    % Get the previous estimate
    Z_prev = Z_hat(:,i-1);

    % Calculate the new estimate
    Z_new = B * AT * Z_prev; % Big A matrix

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

plotTestAccuracy(Z_hat, trans, nu);

%% Define the Mode Process Vector

% Initialise an empty 2 x n_samples matrix to store the mode process
% variable
M_hat = zeros(2, n_samples);

% Use indexing to calculate M_k = M(Z_k)
M_hat(1,:) = Z_hat(1,:);
M_hat(2,:) = 1 - Z_hat(1,:);

%% Cost Function Stopping Time

% Define the penalty each time step
c = 0.001;

% Start by evaluating the summation term for each time step (tau)
% The stopping time must be bounded by n_samples as there is no possibility
% the stopping time can exceed n_samples
J = zeros(1,n_samples);

% Calculate the cost function at each stopping time by looping aruond the
% maximum possible stopping time values 1:n_samples
for i = [1:n_samples]
    J(i) = c*sum(M_hat(2,1:i)) + M_hat(1,i);
end

% Determine the expected value at each tau step
%J = J .* [1:n_samples];

% Determine the value at which J is minimised
[~,tau_J] = min(J);

%% Infimum Bound Stopping Time

% Define a probability threshold to test for
h = 0.99;

% Form a set of k values
k_h = [1:n_samples];

% Get the mode statistic for when the system is in post-change
M_h = M_hat(2,:);
% Index the set such that it contains entries of M^2 > h
k_h = k_h(M_h > h);

% Determine the lower bound of the M_h set to determine the first location
% of when the test statistic exceeds the probability threshold
tau = min(k_h);

% Cleanup
clearvars k_h M_h

%% Plot the stopping results
plotStoppingResults(n_samples,nu,tau,M_hat,h);

%% Calculate performance parameters

% Average Detection Delay
ADD = max(0,tau - nu);

% Probability of False Alarm
PFA = 1 - M_hat(2,tau);

%% Cleanup
clearvars i j