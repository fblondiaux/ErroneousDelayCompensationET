% This function plots data based on the provided inputs.
%
% Parameters:
%   - forces: array of forces
%   - tv: time vector
%   - Y1: data array 1
%   - Y2: data array 2
%   - index: index value
%   - yLabel: label for the y-axis
%   - titleText: title for the plot
%   - thickLine: line width
%   - axPos: position of the subplot
%   - subplotParams: subplot parameters
%   - col_c: color for data array 1 (optional)
%   - col_p1: color for data array 2 (optional)
%   - Y3: data array 3 (optional)
% Returns:
%   - ax: subplot handle
function [ax] = plotData(forces, tv, Y1, Y2, index, yLabel, titleText, thickLine, axPos, subplotParams, col_c, col_p1, col_p2, Y3)
    % Create subplot
    ax = subplot(subplotParams{:}, 'Units', 'centimeters');
    ax.Position = axPos; % define your position
    hold on;

    %Plot the data
    for f = forces
        if nargin == 10
            [col_c, col_p1] = assignColors(f);
        end
        
        plot(tv, squeeze(mean(Y1(f, :, index, :), 2)) * 180 / pi, 'Color', col_c, 'LineWidth', thickLine); %From rad to deg
        plot(tv, squeeze(mean(Y2(f, :, index, :), 2)) * 180 / pi, 'Color', col_p1, 'LineWidth', thickLine);

        if nargin == 14
            plot(tv, squeeze(mean(Y3(f, :, index, :), 2)) * 180 / pi, 'Color', col_p2, 'LineWidth', thickLine);
        end
    end

    xline(0);
    xlabel('Time (ms)');
    ylabel(yLabel);
    title(titleText);
    xlim([tv(1) tv(end)])
end
