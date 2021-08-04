function plotObservationData(n_sensors, ...
    trans, y, nu, mean_unaffected)
% Plot the observation data vs. samples with the number of sensors,
% transition points of the state vector, observation data and changepoint

% Determine the number of samples
n_samples = size(y,2);

figure
y_lim = [-1.2*max(y,[],'all') 1.2*max(y,[],'all')];
for i = [1:n_sensors]
    % Plot the observation vectors and their associated transition points
    subplot(n_sensors,1,i)
    hold on

    % Use the transitions for sensor i to plot where they are occurring on
    % the observation plots
    % Plot the regions of when the sensor nodes are affected by a different
    % statistical distribution as a faded red rectangle
    
    % Get the transitions related to sensor i
    node_ind = trans(:,1) == i+1;
    node_trans = trans(node_ind,:);
    
    % Loop around the transitions for sensor i
    for j = [1:size(node_trans,1)]
        % Fetch the required indexes
        cur_start = node_trans(j,2);
        cur_stop = node_trans(j,3);
        
       % Plot a rectangle that overlays onto the transition points
       rectangle('Position',[cur_start y_lim(1) ...
             cur_stop-cur_start y_lim(2)-y_lim(1)], ...
            'FaceColor',[1 0 0 0.3])
    end
    
    plot([1:n_samples], y(i,:),'b') % Observation vector plot
    xline(nu,'g-') % System changepoint identifier
    yline(mean_unaffected(i),'y--') % Expected mean value
    hold off

    title(['Gaussian Observation y vs. Samples k of sensor ' num2str(i)])
    xlabel('Sample k')
    ylabel('Observation y')
    ylim(y_lim)
end
end

