function [transitions] = getTransitionIndices(e_trans,n_sample_change)
% Fetch the transitions of the DTMC defined in e_trans using a maximum
% sampling count of n_sample_change

% Initialise the transitions variable
% Stored as follows:
% SENSOR NODE || TRANSITION START || TRANSITION STOP
transitions = zeros(n_sample_change,size(e_trans,1));

% Begin looping around the transition indicator matrix
i = 1;
insert = 1;
while i <= n_sample_change
    % Get the current indicator column vector from thet matrix
    e_cur = e_trans(:,i);
    index_cur = max(e_cur);

    % Check whether the current index holds a non-zero transition element
    if index_cur > 0
       affected_start = index_cur;

       % Check whether the next element holds another non-zero element
       % would signify a transition lasting for k=1 sample
       if i < n_sample_change && max(e_trans(:,i+1)) > 0
           affected_stop = affected_start;
           i = i + 1; % Increment by one to move to the index next to this singular element
       else
           % If there is a transition greater than 1 sample in length, keep
           % looping until the next changepoint is found
           j = i + 1;
           while j < n_sample_change && max(e_trans(:,j)) == 0
               if j >= n_sample_change
                   j = n_sample_change; % Check when bounds of array have been exceeded
                   break;
               end
               j = j + 1;
           end

           % Set the endpoint of the transition
           affected_stop = max(e_trans(:,j));

           % Update i to start at the increment after the new transition
           i = j + 1;
       end

       % Check when bounds of array have been exceeded
        if j >= n_sample_change
            i = n_sample_change;
            affected_stop = n_sample_change;
        end

       % Get the row index of the current transition values
       [~,row_cur] = max(e_trans(:,j));

       % Update the transitions matrix
       transitions(insert,:) = [row_cur affected_start affected_stop];
       insert = insert + 1;
       
    else
        % Increment if nothing found
        i = i + 1;
    end
end

% Trim the rest of the transitions matrix away
i = 1;
while i <= n_sample_change && transitions(i,1) > 0
    i = i + 1;
end
transitions(i:n_sample_change,:) = [];
end % Function End

