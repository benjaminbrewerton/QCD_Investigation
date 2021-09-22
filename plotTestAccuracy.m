function plotTestAccuracy(Z_hat,trans, nu)
% Plot the accuracy of the test statistic versus the samples, with the
% affected areas of the observations shown.
% Takes the test statistic, transition points and the changepoint sample
% as inputs.

% Define the number of states and samples
n_states = size(Z_hat,1);
n_samples = size(Z_hat,2);
if length(nu) > 1
    n_sensors = size(Z_hat,1) / 2;
end

figure
y_lim = [-0.25 1.25];
for i = [1:n_states]
    % Plot the test statistics and their associated transition points
    subplot(n_states,1,i)
    hold on

    % Use the transitions for sensor i to plot where they are occurring on
    % the test statistic plots
    % Plot the regions of when the sensor nodes are affected by a different
    % statistical distribution as a faded red rectangle
    
    if i ~= 1
        % Get the transitions related to sensor i
        node_ind = trans(:,1) == i;
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
    end

    
    plot([1:n_samples], Z_hat(i,:),'b') % Test statistic plot
    
    % Loop around the changepoint vector and plot each
    for j = [1:length(nu)]
        if nu(j) < 0
            xline(abs(nu(j)),'g-') % System changepoint identifier
        else
            xline(abs(nu(j)),'r-') % System changepoint identifier
        end
    end

    set(gca, 'color', [0 0.07 0.1 0.2])
    title(['Test Statistic $$Z_k^' num2str(i) '$$ vs. Samples k'],'Interpreter','Latex')
    xlabel('Sample k','Interpreter','Latex')
    ylabel(['$$Z_k^' num2str(i) '$$'],'Interpreter','Latex')
    xlim([0 n_samples])
    ylim(y_lim) % Leave some space in between the top and bottom y-lims
    
    hold off
end
end

