
% plotErrorBar - Plots error bars for multiple datasets.
%
% Parameters:
%   - forces: vector of forces.
%   - Y1: matrix of data for the first dataset.
%   - col1: color for the first dataset.
%   - Y2: matrix of data for the second dataset.
%   - col2: color for the second dataset.
%   - thickLine: line width for the error bars.
%   - axPos: position of the subplot.
%   - subplotParams: parameters for the subplot.
%   - yLabel: label for the y-axis.
%   - titleText: title for the plot.
%   - indexRange: range of indices for the data.
%   - axislimits: limits for the plot axes.
%   - Y3: matrix of data for the third dataset (optional).
%   - col3: color for the third dataset (optional).
% Returns:
%   - ax: subplot handle
function [ax] = plotErrorBar(forces, Y1, col1, Y2, col2, thickLine, axPos, subplotParams, yLabel, titleText, indexRange, axislimits, Y3, col3)
    % Create subplot
    ax = subplot(subplotParams{:}, 'Units', 'centimeters');
    ax.Position = axPos; % define your position
    hold on;

    %Plot the data
    errorbar(forces, squeeze(mean(mean(Y1(:, :, 1, indexRange), 4), 2))', [0 0 0], '-o', 'Color', col1, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', col1, 'MarkerSize', 5);
    errorbar(forces, squeeze(mean(mean(Y2(:, :, 1, indexRange), 4), 2))', [0 0 0], '-o', 'Color', col2, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', col2, 'MarkerSize', 5);

    if nargin == 14
        errorbar(forces, squeeze(mean(mean(Y3(:, :, 1, indexRange), 4), 2))', [0 0 0], '-o', 'Color', col3, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', col3, 'MarkerSize', 5);
    end

    xlabel('Perturbation (Nm)');
    ylabel(yLabel);
    title(titleText);
    axis(axislimits)
end
