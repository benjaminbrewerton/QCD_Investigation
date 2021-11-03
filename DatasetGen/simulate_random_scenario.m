function [X, y, nu] = simulate_random_scenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected,nu)

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

% Check for a deterministic changepoint setting
if ~exist('nu','var')

% Simulate the entire scenarios markov chain state transitions
% Assume that the when the space transitions from alpha to beta that the
% simulation begins at a random state with equal probability of initial
% state
% Initialise the sequence to begin at node 1
X = simulate(A, n_samples - 1,'X0',[1 zeros(1,n_sensors)]);

% Determine nu (changepoint) as the first time the sequence escapes from
% DTMC alpha and into DTMC beta
nu = length(X(X == 1)) + 1; % + 1 for being inclusive of the transition sample

while nu >= n_samples
    disp(['Changepoint was never reached after ' num2str(n_samples) ...
       ' samples. Trying again.']);
    % Initialise the sequence to begin at node 1
    X = simulate(A, n_samples - 1,'X0',[1 zeros(1,n_sensors)]);
    nu = length(X(X == 1)) + 1; % + 1 for being inclusive of the transition sample
end

else % Deterministic case
    % Initialise the state sequence with ones until the changepoint to indicate
    % the system is in the pre-change state
    X = ones(1, n_samples);

    if nu ~= 0
        % Simulate the rest of the markov chain in the post-change state space
        % using the post-change transition matrices
        X(nu:end) = simulate(A_beta, n_samples - nu) + 1;
    end
end
%% Generate randomly distributed values for each sensing node

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors,n_samples);

% Populate the y matrix with samples taken from a Gaussian distribution
% with mean and variances per state as defined in means and vars
for i = [1:n_samples]
    % Check whether we are in pre or post-change
    if(i < nu || nu == 0)
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

end