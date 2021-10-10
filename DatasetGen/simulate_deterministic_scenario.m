function [X, y, nu] = simulate_deterministic_scenario(mean_unaffected,var_unaffected, ...
    mean_affected,var_affected, nu)

%% Begin definition of network variables

% Number of sensors for post-change modelling
n_sensors = 3;

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
if ~exist('nu','var')
    % Generate a random changepoint between 0 and n_samples / 2
    nu = randi([1 n_samples-1]);
end

% Check if nu is invalid
while nu >= n_samples
    disp("Regenerating changepoint of deterministic scenario")
    nu = randi([1 n_samples]);
end

% Initialise the state sequence with ones until the changepoint to indicate
% the system is in the pre-change state
X = ones(1, n_samples);

if nu ~= 0
    % Simulate the rest of the markov chain in the post-change state space
    % using the post-change transition matrices
    X(nu:end) = simulate(A_beta, n_samples - nu) + 1;
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