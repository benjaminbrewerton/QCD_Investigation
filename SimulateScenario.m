%% Environment Clearing
clc
close all
clear vars

%% Begin definition of network variables

% Number of sensors
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
A_nu = ones(1,n_sensors) ./ n_sensors;

% Plot the Beta Sensor space transition probabilities
figure
graphplot(A_beta,'ColorEdges',true,'LabelEdges',true)
title('Post-Change Transition Probabilities')

%% Define the probabilities of specific events

% There is no pre-change state dependence, such that rho is constant for
% all states in the state spaces

% Rho is the probability that an object will enter the network
rho = 5e-4;

% Define the geometric prior to represent the probability of the state not changing
pi = (1-rho).^((1:n_samples) - 1) * rho;

% The probability that the network will tranisition from state alpha to
% beta
P_v = rho * A_nu;

%% Generate the state randoms and changepoint

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

% Determine nu (changepoint) as the first time the sequence escapes from
% DTMC alpha and into DTMC beta
nu = length(X_sys(X_sys == 1)) + 1; % + 1 for being inclusive of the transition sample
% Print the result
disp(['It took ' num2str(nu) ' iterations to transition to space beta']);

% Post-change sample count
n_sample_change = n_samples - nu;

% Assume that the when the space transitions from alpha to beta that the
% simulation begins at a random state with equal probability of initial
% state
X = simulate(A_beta, n_sample_change - 1);

%% Determine the transition points

% Define the e matrix as being 1 when the sensor node at row j and 0
% everywhere else at time k
e = zeros(n_sensors,n_sample_change);

% Use a for loop to iterate around the sample vector and update e to
% determine what sensor's transitions happen at what sample time k
for i = [1:n_sample_change]
    e(X(i),i) = 1;
end

% Determine the sample index of each transition so that it can be
% visualised later
e_trans = zeros(n_sensors,n_sample_change);

% Loop to isolate each indicator vector index
last_sensor = 0;
for i = [1:n_sample_change]
    % Search for which row (sensor) index that the indicator vector shows
    [~,e_node] = min(abs(e(:,i) - 1));
    if last_sensor ~= e_node
        e_trans(e_node,i) = i; % Set the transition to this particular k index
        last_sensor = e_node; % Update the last known sensor change
    end
end

% Cleanup
clearvars e_node last_sensor

%% Define the Normal Distributions for each sensing node

% Assume the sensors observation measurements are i.i.d.

% Each sensor's default distribution without being affected by a target
% will be a standard normal distribution with sigma = 1 and mean = 0
dist_mean = 0;
dist_dev = 1;

% The modifier for how much the variance should be increased by when the
% sensor node is affected by the target object
stat_modifier = 2;

% Each distribution parameter can be accessed with the node's index
variances = 0.25 + 0.5*rand(1,n_sensors);

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors,n_sample_change);

for i = [1:n_sensors]
    % Define the non-affected sensor's distribution
    y(i,:) = random(makedist('normal',dist_mean,dist_dev),1,n_sample_change);
    
    % Initialise the affected distribution for each sensor node
    current_dist = makedist('normal',2,variances(i)*stat_modifier);
    
    % Get the current e_trans vector
    e_cur = e_trans(i,:);
    % Trim the vector such that only non-zero elements remain
    e_cur(e_cur == 0) = [];
    
    % Loop around e_cur and determine where to insert the modified samples
    for j = 1:2:length(e_cur)
       % Determine the stop and start points of the current transition
       affected_start = e_cur(j);
       if j == length(e_cur)
           affected_stop = n_sample_change;
       else
           affected_stop = e_cur(j + 1) - 1;
       end
       
       % Update the distribution between the affected points
       y(i,affected_start:affected_stop) = random(current_dist,1, ...
            affected_stop - affected_start + 1);
    end
end

% Plot the histogram of one of the sensor's generated samples
sensor_plot = randi([1 n_sensors]);
figure
histfit(y(sensor_plot,:))
title(['Histogram of sensor ' num2str(sensor_plot)])
ylabel('Sample Count')
xlabel('Value Magnitude')

% Also plot the generated samples vs. sample iteration
figure
y_lim = [-4 4];
for i = [1:n_sensors]
    % Get the current transition row
    e_cur = e_trans(i,:);
    % Trim the transition row to only include non-zero values
    e_cur(e_cur == 0) = [];
    
    % Plot the observation vectors and their associated transition points
    subplot(n_sensors,1,i)
    hold on
    
    % Plot the regions of when the sensor nodes are affected by a different
    % statistical distribution as a faded red rectangle
    for j = [1:2:length(e_cur)]
       % Determine the stop and start points of the current transition
       affected_start = e_cur(j);
       if j == length(e_cur)
           affected_stop = n_sample_change;
       else
           affected_stop = e_cur(j + 1);
       end
       
       % Plot a rectangle that overlays onto the transition points
       rectangle('Position',[affected_start y_lim(1) ...
             affected_stop-affected_start y_lim(2)-y_lim(1)], ...
            'FaceColor',[1 0 0 0.35])
    end
    
    plot([1:n_sample_change], y(i,:),'b') % Observation vector plot
    plot(e_cur(e_cur > 0),0,'g.') % Transition Points
    hold off

    title(['Gaussian Observation y vs. Samples k of sensor ' num2str(i)])
    xlabel('Sample k')
    ylabel('Observation y')
    xlim([0 n_sample_change])
    ylim(y_lim)
end

% Clean up the workspace
clearvars sensor_plot current_dist %e_cur

%% Determine the probability of each observation

% Calculate the z-values of all the probabilities in the observation table
Z = zscore(y,1,2);

% Calculate the probability that each z-value has of being less than or
% equal to a randomly generated variable X
P_z = normcdf(Z);
