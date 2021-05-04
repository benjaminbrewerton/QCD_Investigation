function [colour_legend] = generateColourLegend(cols,text)
% Takes a colour input matrix (cols) with each row being a set of RGB values and
% assigns them to a textual legend reference in a cell array format

% Check for mismatching in the legend entries
if size(cols,1) ~= length(text)
    disp("Legend Size Mismatch");
else
    % Initialise the cell array which will contain the legend entries
    colour_legend = [];
    
    % Loop over all the colour entries in the cols matrix
    for i = [1:length(text)]
        colour_legend = [colour_legend strcat('{\color[rgb]{', ...
            num2str(cols(i,1)), " ", num2str(cols(i,2)), " ", ...
            num2str(cols(i,3)), " 0.3", "} ■ } ", text(i))];
        
        %{'{\color{red} o } Red', '{\color{blue} o } Blue', ...
    %'{\color{black} o } Black'}
    end
end
end
