function [ADD,PFA,tau,S] = ISD_Deterministic(mean_mu,var_mu, ...
    mean_nu,var_nu, y, lambda, h)
%% Begin definition of network variables

% Number of pre-change states
n_states_mu = 1;

% Number of states in the post-change state
n_states_nu = 3;

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

%% Infimum Bound Stopping Time

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
end

% Cleanup
clearvars M_stat

%% Calculate performance parameters

% Average Detection Delay
ADD = mean(abs(abs(lambda) - tau));

% Probability of False Alarm
PFA = 0;

end