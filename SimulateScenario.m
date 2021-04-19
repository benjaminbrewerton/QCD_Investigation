%% Environment Clearing
clc
close all
clear vars

%% Begin definition of network variables

% Number of sensors
n_sensors = 3;

% Number of random samples
n_samples = 1e5;

% Generate the primary changepoint of the system which will transition the
% state from one to another
% Generate a random variable between 25% and 50% of the sample range
nu_initial = randi([n_samples/4 n_samples/2]);

% Generate a random changepoint (nu) for each sensing node that falls a
% minimum of 1000 samples from the primary system changepoint
nu_sensors = zeros(1,n_sensors);
for i = [1:n_sensors]
    nu_sensors(i) = randi([nu_initial+1e3 3/4*n_samples]);
end

% Override one of the sensors to be affected by the primary changepoint as
% would be expected by an object entering a sensor network and one of the
% sensors being affected by this same statistical change
nu_sensors(randi([1 n_sensors])) = nu_initial;

% Generate the sensor's change order to know which sensor's changepoint is
% first affected, then the subsequent changepoints after that
[~,nu_order] = sort(nu_sensors);

% Generate a nu matrix which defines the changepoint (column 1) and ends of
% statistical change (column 2)
nu = zeros(n_sensors,2);

%% Define Sensor States

% Generate an empty array for each sensing node of pre-change and post-change 
% distributions
S_alpha = [];
S_beta = zeros(n_sensors);

% Since there is only one sensing node modelled in the pre-change
% distribution, the indicator vector can be easily inserted
S_alpha = [S_alpha generateIndicatorVector(1,1)];

% There are n_sensors nodes being modelled in post-change distribution beta.
% To initialise these nodes, a for loop will be used for each sensor
for i = [1:n_sensors]
   S_beta(:,i) = generateIndicatorVector(i,n_sensors);
end

%% Define Transition Matrices for states

% Form the Discrete Time Markov Chains for each sensor space

% Since in state alpha, a singular non-sensed state is modelled:
A_alpha = dtmc([1]);

% Define a transition modifier which models what the probability is of a
% detected object moving outside of the current state
P_move = 0.02;

% In state beta, there are n_sensors nodes, so to create a random
% transition matrix
% Uncomment the line below for a fully randomised state transition matrix
%A_beta = mcmix(n_sensors);

% This beta DTMC is formed from the P_move variable on the diagonal and all
% other values are P_move/2
trans_beta = diag(repelem(1 - P_move, n_sensors));
trans_beta(trans_beta == 0) = P_move / 2;
A_beta = dtmc(trans_beta);

% For transition state nu, assume that the probability for the state to
% tranisition to the states in space beta is equal
A_nu = ones(1,n_sensors) ./ n_sensors;

% Define the M_X vector as being 1 when the changepoint has been reached and
% 0 when it has not yet been reached for each sensor sampling
M = zeros(n_sensors,n_samples);

% Plot the Alpha Sensor space transition probabilities
figure
graphplot(A_alpha,'LabelEdges',true)
title('Pre-Change Transition Probability')
% Plot the Beta Sensor space transition probabilities
figure
graphplot(A_beta,'ColorEdges',true,'LabelEdges',true)
title('Post-Change Transition Probabilities')

%% Define the Normal Distributions for each sensing node

% Assume the sensors observation measurements are i.i.d.

% Each sensor's default distribution without being affected by a target
% will be a standard normal distribution with sigma = 1 and mean = 0
dist_mean_unaffected = 0;
dist_dev_unaffected = 1;

% Each distribution parameter can be accessed with the node's index

% Define the means of each sensor, here we just use 0 for each sensor
dist_mean = repelem(0,n_sensors);

% Define the standard deviations of each sensor using a randomly generated
% distribution which are 1 <= sigma <= 2
dist_devs = 1 + rand(1,n_sensors);

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors,n_samples);

for i = nu_order([1:n_sensors])
    % Get the current iteration of nu_order
    [~,j] = min(abs(nu_order - i));
    
    % Define the non-affected sensor's distribution
    y(i,:) = random(makedist('normal',dist_mean_unaffected,dist_dev_unaffected),1,n_samples);
    
    % Initialise the affected distribution for each sensor node
    %current_dist = makedist('normal',dist_mean(i),dist_devs(i));
    stat_modifier = 1.5;
    current_dist = makedist('normal',dist_mean_unaffected,dist_dev_unaffected*stat_modifier);
    
    % Define where in the distributions to insert the affected node
    % distribution samples
    y_a_start = nu_sensors(i);
    y_a_stop = -1;
    
    % Use an if statement to check whether there is another changepoint
    % beyond the current one in the system
    if j < n_sensors
        y_a_stop = nu_sensors(nu_order(j+1));
    else
        y_a_stop = n_samples + 1; % + 1 to account for the sample overlap prevention on L~126
    end
    
    % Update the nu matrix to use the new changepoint stop points
    nu(j,:) = [y_a_start y_a_stop];
    
    % Update M to reflect the changepoint for each node
    M(j,y_a_start:y_a_stop) = 1;
    
    disp(['sensor: ' num2str(i) ' at index: ' num2str(j) ', change: ' num2str(y_a_start) ', stop: ' num2str(y_a_stop-1)]);
    
    % Pull n_samples worth of randomised data from the established
    % distribution for when the sensor is affected by a target
    y(j,y_a_start : y_a_stop - 1) = random(current_dist,1,y_a_stop - y_a_start);
end

% Define another observation variable y_nu which contain the observations
% taken after the changepoint has occurred
y_nu = y(:,nu:n_samples);

% Plot the histogram of one of the sensor's generated samples
sensor_plot = randi([1 n_sensors]);
figure
histfit(y(sensor_plot,:))
title(['Histogram of sensor ' num2str(sensor_plot)])
ylabel('Sample Count')
xlabel('Value Magnitude')

% Also plot the generated samples vs. sample iteration
figure
for i = [1:n_sensors]
    subplot(n_sensors,1,i)
    hold on
    plot([1:n_samples], y(i,:)) % Observation vector plot
    plot(nu(i,1),0,'r.') % Changepoint
    plot(nu(i,2),0,'g.') % End Changepoint
    hold off
    title(['Gaussian Observation y vs. Samples k of sensor ' num2str(i)])
    xlabel('Sample k')
    ylabel('Observation y')
    xlim([0 n_samples])
end

% Clean up the workspace
clearvars sensor_plot current_dist

%% Determine the probability of each observation

% Calculate the z-values of all the probabilities in the observation table
Z = zscore(y,1,2);

% Calculate the probability that each z-value has of being less than or
% equal to a randomly generated variable X
P_z = normcdf(Z);

%% Define the probabilities of specific events

% There is no pre-change state dependence, such that rho is constant for
% all states in the state spaces

% Rho is the probability that an object will enter the network
rho = 5e-4;

% Define the geometric prior to represent the probability of the state not changing
pi = (1-rho).^((1:n_samples) - 1) * rho;

% The probability that the network will tranisition from state alpha to
% beta
P_v = rho *A_nu;

%% Generate the state randoms

% Since before the changepoint, the state variables in space alpha are not
% relevant to the problem and do not have to be simulated.

% Calculate A to determine the state variables X
% Use the definition where rho is constant for all states
A = dtmc([(1-rho)*A_alpha.P P_v ; zeros(n_sensors,1) A_beta.P], ...
    'StateNames',["Alpha 1" "Beta 1" "Beta 2" "Beta 3"]);

% Display the transition probabilities of the overall system
figure
graphplot(A,'ColorEdges',true)
title('Holistic System Transition Probabilities')

% Simulate the entire scenarios markov chain state transitions
X_sys = simulate(A, n_samples - 1,'X0',[1 zeros(1,n_sensors)]);
samples_escape = length(X_sys(X_sys == 1)) + 1; % + 1 for being inclusive of the transition sample
% Print the result
disp(['It took ' num2str(samples_escape) ' iterations to transition to space beta']);

% Define the M_X vector as being 1 when the changepoint has been reached and
% 0 when it has not yet been reached
M_X = zeros(1,n_samples);
M_X(X_sys > 1) = 1;

% Assume that the when the space transitions from alpha to beta that the
% simulation begins at state 1
X = simulate(A_beta, n_samples - nu_initial);

% Calculate the mode process transition dtmc
% A_M = dtmc([1-rho rho ; 0 1]);
% Simulate A_M
% M = simulate(A_M, n_samples);
