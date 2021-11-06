function plotStoppingResults(nu, tau, M_hat, h)
% Plot the stopping results of the QCD scheme with arguments: number of
% samples, changepoint sample index, stopping time, posterior probability
% check (M) and detection threshold

% Declare the number of samples in the scenario
n_samples = size(M_hat,2);

% Take the absolute of the changepoints
nu_abs = abs(nu);

if length(nu) > 1 || length(tau) > 1 % Multiple changepoints
    % Define the titles
    titles = {'Pre-change Mode Process $$\hat{M}_k^\mu$$ vs. Samples k', ...
        'Post-change Mode Process $$\hat{M}_k^\nu$$ vs. Samples k'};
    legends = {'$$\hat{M}_k^\mu$$ -- $$P(X\in S_{\mu})$$', ...
        '$$\hat{M}_k^\nu$$ -- $$P(X\in S_{\nu})$$'};
    ytitles = {'\mu','\nu'};
    
    figure
    for j = [1:size(M_hat,1)]
        subplot(2,1,j)
        hold on
        
        % Declare an index vector for only displaying the correct
        % changepoints on the likelihood sets
        if j == 1
           trim = nu < 0; 
        else
           trim = nu > 0;
        end
        
        % Trim the stopping times that are null away
        tau_trim = tau >= 0;
        
        plot([1:n_samples],M_hat(j,:)) % Plot likelihood of system in post-change
        
        plot(nu_abs(trim),M_hat(j,nu_abs(trim)),'go') % Plot the change-points

        plot(tau(trim & tau_trim),M_hat(j,tau(trim & tau_trim)),'ro') % Plot the stopping times
        
        yline(h,'m--') % Plot the detection threshold

        set(gca, 'color', [0 0.07 0.1 0.2])
        title(char(titles(j)),'Interpreter','Latex')
        xlabel('Sample k','Interpreter','Latex')
        ylabel(['$$M_k^' ytitles{j} '$$'],'Interpreter','Latex')
        leg = legend(char(legends(j)),...
            '$$\nu$$ -- Changepoints', '$$\tau$$ -- Stopping Times', ...
            '$$h$$ -- Threshold');
        leg.Interpreter = 'Latex';
        leg.Color = 'w';
        xlim([0 n_samples])
        ylim([-0.1 1.1])
    end
    
else % Singular changepoint
    figure
    hold on

    plot([1:n_samples],M_hat) % Plot likelihood of system in post-change
    xline(nu,'g--') % Plot the change-points
    plot(tau,M_hat(tau),'ro') % Plot the stopping time
    yline(h,'m--') % Plot the detection threshold

    set(gca, 'color', [0 0.07 0.1 0.2])
    title('Post-change Mode Process $$\hat{M}_k^2$$ vs. Samples k','Interpreter','Latex')
    xlabel('Sample k','Interpreter','Latex')
    ylabel('$$\hat{M}_k$$','Interpreter','Latex')
    leg = legend('$$\hat{M}_k$$ -- $$P(Z_k\in S_{\nu})$$',...
        '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time', ...
        '$$h$$ -- Threshold');
    leg.Interpreter = 'Latex';
    leg.Color = 'w';
    xlim([0 n_samples])
    ylim([-0.1 1.1])
end
end

