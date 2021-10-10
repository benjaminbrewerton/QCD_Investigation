function [transitions] = getTransitionIndices(X_seq)
% Fetch the transitions of the DTMC defined in e_trans with optional change
% limit defined by delta_lim.

% Get the samples used in e_trans
n_samples = size(X_seq,1);

% Initialise the transitions variable
% Stored as follows:
% SENSOR NODE || TRANSITION START || TRANSITION STOP
transitions = zeros(n_samples,3);

% Define holder variables for the previous node
prev_node = 1;
prev_trans = 1;
prev_change = 1;

% Begin looping around the transition indicator matrix
for i = [1:n_samples]
    cur_node = X_seq(i);
    if cur_node ~= prev_node
        % Update the current found transition
        transitions(prev_trans,:) = [prev_node prev_change i-1]; 
        % Update the previous node
        prev_node = X_seq(i);
        % Increment the previous transition
        prev_trans = prev_trans + 1;
        % Update the previous changepoint
        prev_change = i;
    end
end

% Check for overflow
if(prev_change+1 >= 10000)
    prev_change = n_samples - 1;
end

if prev_change > n_samples - 1
    % Rare occasion when transition point is at k=9999
    transitions(prev_trans,:) = [X_seq(n_samples) prev_change n_samples];
else
    % Determine the transition point of the remainder sample sequence
    final_change = transitions(prev_trans - 1, 3);
    final_node = X_seq(final_change + 1);

    % Update the final row of the transition matrix to reflect the skipped
    % transition point determined previously
    transitions(prev_trans,:) = [final_node prev_change n_samples];
end

% Trim the first row away as it represents the target when it is still not
% encountered the changepoint
%transitions(1,:) = [];

% Trim the rest of the transitions matrix away
transitions_ind = transitions(:,1) == 0;
transitions(transitions_ind,:) = [];

end % Function End

