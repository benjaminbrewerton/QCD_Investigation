function [labelArray] = generateAxisLabels(str_format,label_count)
% Generate axis labels with a set tick spacing of 1 graph unit. Outputs a
% cell array with a space between the top and bottom of the targeted axis.
% Specify str_format with an '~' where the numbering should go in the
% string

% Initialise the output label array
labelArray = cell(1,label_count + 2); % +2 for axis padding

% Loop around the numbering array
for i = [1:label_count + 2]
    if i == 1 || i == label_count+2 % Check for special labelling
        labelArray(i) = {''};
    else
        labelArray(i) = {strrep(str_format,'~',num2str(i-1))};
    end
end

end

