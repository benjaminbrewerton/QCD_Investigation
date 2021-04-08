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
nu = randi([n_samples/4 n_samples/2]);

% Generate a random changepoint (nu) for each sensing node that falls a
% minimum of 1000 samples from the primary system changepoint
nu_sensors = zeros(1,n_sensors);
for i = [1:n_sensors]
    nu_sensors(i) = randi([nu+1e3 3/4*n_samples]);
end

% Override one of the sensors to be affected by the primary changepoint as
% would be expected by an object entering a sensor network and one of the
% sensors being affected by this same statistical change
nu_sensors(randi([1 n_sensors])) = nu;


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

%% Define Transition Matrices for state Alpha and Beta

% Form the Discrete Time Markov Chains for each sensor space
% Since in state alpha, only a singular sensor is modelled:
A_alpha = dtmc([1]);

% In state beta, there are n_sensors nodes, so to create a random
% transition matrix
A_beta = mcmix(n_sensors);

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

% Each distribution parameter can be accessed with the node's index

% Define the means of each sensor, here we just use 0 for each sensor
dist_mean = repelem(0,n_sensors);

% Define the standard deviations of each sensor using a randomly generated
% distribution
dist_devs = rand(1,n_sensors);

% Generate randomised observation data derived from the normal
% distributions for each sensor node. Store the data as a matrix with each
% row being a sensor's observation data and the column's being the samples
y = zeros(n_sensors,n_samples);
for i = [1:n_sensors]
    % Initialise the distribution for each sensor node
    current_dist = makedist('normal',dist_mean(i),dist_devs(i));
    
    % Pull n_samples worth of randomised data from the established
    % distribution
    y(i,:) = random(current_dist,1,n_samples);
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
plot([1:n_samples], y(sensor_plot,:))
title(['Generated Values vs. Sample Count of sensor ' num2str(sensor_plot)])
xlabel('Sample Count')
ylabel('Value Magnitude')

clearvars sensor_plot

%% Determine the probability of each observation

% Calculate the z-values of all the probabilities in the observation table
Z = zscore(y,1,2);

% Calculate the probability that each z-value has of being less than or
% equal to a randomly generated variable X
P_z = normcdf(Z);

%% Generate the States of the Markov Chain
