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

% Since in state alpha, only a singular sensor is modelled:
A_alpha = 1;

% In state beta, there are n_sensors nodes, so to create a random
% transition matrix
A_beta = mcmix(n_sensors);