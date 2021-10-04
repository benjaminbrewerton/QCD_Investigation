function [X, y, lambda] = simulate_ISD_scenario(mean_mu,var_mu, ...
    mean_nu,var_nu, changepoints, lambda)

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

%% Generate the state randoms and changepoint

% Check whether the changepoint is already supplied
if ~exist('nu','var')
    if changepoints < 0 % Random dataset
        % Bayesian Probability Space Assumptions

        % These variables are used for Bayesian changepoint modelling.

        % Rho is the probability that an object will transition between state
        % spaces.
        rho = 2.5e-3;

        % The probability that the network will tranisition from state alpha to
        % beta
        A_mu = rho .* pi_mu;
        A_nu = rho .* pi_nu;

        % Generate the state randoms and changepoint

        % Define the A from literature
        A_lit = [(1-rho).*A_alpha.P (1/n_states_nu).*repelem(A_mu,n_states_mu,n_states_nu) ; ...
            repelem(A_nu.',n_states_mu,1) (1-rho .* (1/n_states_nu)).*A_beta.P];

        % Determine the large A matrix which governs 
        A = dtmc(A_lit.');

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
    elseif changepoints > 0 % Deterministic dataset
        % Generate the checker pattern
        checkers = ones(1,changepoints);
        for i = [1:changepoints]
            if mod(i,2) == 0
                checkers(i) = -1;
            end
        end
        
        % Generate the changepoint set
        lambda = sort(randi([1 n_samples],1,changepoints)) .* checkers;
        
        % Ensure that the maximum changepoint value is n_samples-1
        while max(abs(lambda)) >= n_samples || length(unique(abs(lambda))) ~= changepoints
            disp("Regenerating a changepoint...")
            lambda = sort(randi([1 n_samples],1,changepoints)) .* checkers;
        end
        
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
                if(lambda_next-lambda_cur == 0)
                    disp(lambda)
                end
                X(lambda_cur:lambda_next) = simulate(A_beta, lambda_next-lambda_cur) + 1;
            end
        end
    else % Null changepoint scenario
        
        % Initialise the state sequence with ones until the changepoint to indicate
        % the system is in the pre-change state
        X = ones(1, n_samples);
        % Nullify the changepoint vector
        lambda = [];
    end
else
    % Handle if the changepoint set is already specified 
   
    % Deterministic if anything else

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
            X(lambda_cur:lambda_next-1) = simulate(A_beta, lambda_next-lambda_cur-1) + 1;
        end
    end
end

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

end