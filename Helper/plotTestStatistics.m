function plotTestStatistics(Z_hat,trans)
% Plot the test statistic for the posterior checks using colour coded plots
% using colours defined in this function. This function assumes that only
% three sensors are present in the post-change event.
% Takes the posterior check Z_hat and transition points as inputs.

% Define the number of states and samples
n_states = size(Z_hat,1);
n_samples = size(Z_hat,2);


% Generate unique colour schemes for the sensors
%colours = rand(n_sensors, 3); % 1 random for each RGB value
colours = [1 0 1 ; 1 0 0 ; 0.08 0.5 0 ; 0 0 1];

% Y limits
y_lim = [0 n_states+1];

% ===== Plot of each sensor's Z test statistic =====
figure

for i = [1:n_states]
    subplot(n_states+1,1,i)
    
    plot([1:n_samples], Z_hat(i,:), 'color', colours(i,:))
    
    set(gca, 'color', [0 0.07 0.1 0.2])
    title(['Test statistic $$\hat{Z}_k^' num2str(i) '$$ vs. Samples k'],'Interpreter','Latex')
    ylabel(['$$\hat{Z}_k^' num2str(i) '$$'],'Interpreter','Latex')
    xlim([0 n_samples])
    ylim([-0.25 1.25]) % Leave some space in between the top and bottom y-lims
    
    % Create the legend on the first subplot
    if i == 1
       % Create the legend
        text(0.8*n_samples, 1/8*max(y_lim), ...
        generateColourLegend(colours, ...
        {'e^\alpha_1', 'e^\beta_1', 'e^\beta_2', 'e^\beta_3'}), ...
        'EdgeColor', 'k', 'BackgroundColor', 'w')
    end
    
    % Change the x-axis to be blank
    xticks([0 n_states+1])
    xticklabels(generateAxisLabels(" ",0))
end

subplot(n_states+1,1,n_states+1)
% Loop around each transition point and overlay them with different colours
% depending on what node the transition occurred at
for j = [1:size(trans,1)]
    % Fetch the required indexes
    cur_node = trans(j,1);
    cur_start = trans(j,2);
    cur_stop = trans(j,3);

    % Plot a rectangle that overlays onto the transition points
    rectangle('Position',[cur_start y_lim(1) ...
        cur_stop-cur_start y_lim(2)-y_lim(1)], ...
        'FaceColor',[colours(cur_node,:)], 'EdgeColor',[0 0 0 0])
end

hold off

% Change the y-axis to be blank
yticks([0 n_states+1])
yticklabels(generateAxisLabels(" ",0))

set(gca, 'color', [0 0.07 0.1 0.2])
title('Generated State Sequence $$X_k$$ vs. Samples k','Interpreter','Latex')
xlabel('Sample k')
ylabel('$$X_{k}$$','Interpreter','Latex')
xlim([0 n_samples])
ylim(y_lim) % Leave some space in between the top and bottom y-lims
end

