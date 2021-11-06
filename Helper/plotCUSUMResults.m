function plotCUSUMResults(S,h,nu,tau)
% Plot the CUSUM stopping results

% Declare the number of samples in the scenario
n_samples = size(S,2);

% Take the absolute of the changepoints
nu_abs = abs(nu);

% Trim the stopping times that are null away
tau_trim = tau >= 0;

figure
hold on

for i = 1:size(S,1)
    plot([1:n_samples],S(i,:)) % Plot likelihood of system in post-change
end

if size(S,1) > 1
    for i = 1:length(nu)
        if nu(i) < 0
            plot(nu_abs(i),S(1,nu_abs(i)),'go') % Plot the change-point
            if tau(i) > -1
                plot(tau(i),S(2,tau(i)),'ro') % Plot the stopping times
            end
        else
            plot(nu_abs(i),S(2,nu_abs(i)),'go') % Plot the change-point
            if tau(i) > -1
                plot(tau(i),S(1,tau(i)),'ro') % Plot the stopping times
            end
        end
    end
    
else
    xline(nu_abs,'g--') % Plot the change-point
    plot(tau(tau_trim),S(tau(tau_trim)),'ro') % Plot the stopping times
end

for i = [1:length(h)]
    yline(h(i),'m--') % Plot the detection threshold
end

set(gca, 'color', [0 0.07 0.1 0.2])
% Only set log scale on traditional QCD scenario
if length(nu) <= 1 || length(tau) <= 1
    set(gca, 'YScale', 'log')
    ylim([0 5000])
end
title('CUSUM Test Statistic $$\tilde{Z}_k$$ vs. Samples k','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$\tilde{Z}_k$$','Interpreter','Latex')

if size(S,1) > 1
    leg = legend('$$S_\nu$$ -- CUSUM Statistic',...
    '$$S_\mu$$ -- CUSUM Statistic',...
    '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time');
else
    leg = legend('$$S_\nu$$ -- CUSUM Statistic',...
    '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time', ...
    '$$h_C$$ -- Threshold');
end
leg.Interpreter = 'Latex';
leg.Color = 'w';
xlim([0 n_samples])

end

