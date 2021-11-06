function [ADD,PFA,tau,M_hat] = ISD_Random(mean_mu,var_mu, ...
    mean_nu,var_nu, y, lambda, h)
% Produces the Average Delay Detection and Probability of False Alarm from
% a set of inputted means and variances (unaffected and affected
% respectively in order) for an ISD scenario. 
% This function assumes a setup of 3 sensors in the network.

%% Begin definition of network variables

% Number of pre-change states
n_states_mu = 1;

% Number of states in the post-change state
n_states_nu = 3;

% Total number of states
n_states = n_states_mu + n_states_nu;

% Number of sensors
n_sensors = max(n_states_mu, n_states_nu);

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

%% Define the Mode Process Vector

% Initialise an empty 2 x n_samples matrix to store the mode process
% variable
M_hat = zeros(2, n_samples);

% Use indexing to calculate M_k = M(Z_k)
M_hat(1,:) = Z_hat(1,:);
M_hat(2,:) = 1 - Z_hat(1,:);

if ~exist('h','var')
    ADD = 0;
    PFA = 0;
    tau = 0;
    return
end

%% Infimum Bound Stopping Time

% Get the mode statistic for when the system is in post-change
M_mu = M_hat(1,:);
M_nu = M_hat(2,:);

% Initialise a vector which stores the threshold crossing indices. Store as
% -1 if no changepoint is detected within the common boundary.
tau = ones(1,length(lambda)) .* -1;

% Loop around the changepoints and determine the index of the first sample
% which exceeds the specified thresholds
if exist('lambda','var')
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
end

% Cleanup
clearvars M_stat

%% Calculate performance parameters

DD = tau - abs(lambda);

% Probability of False Alarm
PFA = sum(DD < 0);

% Average Detection Delay
DD(DD < 0) = n_samples;
ADD = mean(DD);

end