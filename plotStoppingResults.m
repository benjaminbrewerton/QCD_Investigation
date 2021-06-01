function plotStoppingResults(n_samples,nu, ...
    tau, M_hat, h)
% Plot the stopping results of the QCD scheme with arguments: number of
% samples, changepoint sample index, stopping time, posterior probability
% check (M) and detection threshold

figure
hold on

plot([1:n_samples],M_hat(2,:)) % Plot likelihood of system in post-change
plot(nu,M_hat(2,nu),'go') % Plot the change-point
plot(tau,M_hat(2,tau),'ro') % Plot the stopping time
yline(h,'m--') % Plot the detection threshold

set(gca, 'color', [0 0.07 0.1 0.2])
title('Post-change Mode Process $$\hat{M}_k^2$$ vs. Samples k','Interpreter','Latex')
xlabel('Sample k','Interpreter','Latex')
ylabel('$$M_k^2$$','Interpreter','Latex')
leg = legend('$$\hat{M}_k^2$$ -- $$P(X\in S_{\beta})$$',...
    '$$\nu$$ -- Changepoint', '$$\tau$$ -- Stopping Time', ...
    '$$h$$ -- Threshold');
leg.Interpreter = 'Latex';
leg.Color = 'w';
xlim([0 n_samples])
ylim([-0.1 1.1])
end

