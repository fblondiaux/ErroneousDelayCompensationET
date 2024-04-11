% This function saves the figure F as an SVG file using Inkscape. 
% Changing font, legend size, lines thickness, ... before saving the figure
%
% Parameters:
%    - F: the figure handle
%    - fullFileName: the full file path and name of the SVG file to be saved
function figForInkscapeSave(F, fullFileName)

    constantsPlots;

    figure(F);
    set(gca, 'linewidth', thinLine)
    set(gca, 'Color', 'none'); % Sets axes background
    savefig(fullFileName)
    set(gca, 'box', 'off')
    set(findall(gcf, '-property', 'FontSize'), 'FontWeight', 'Normal', 'FontSize', 6, 'fontname', 'Arial')

    % Update the font size of the legend
    legendFontSize = 4;
    hLegend = findobj(gcf, 'Type', 'legend');
    set(hLegend, 'FontSize', legendFontSize);

    % Update the font size of the title
    set(F, 'PaperPositionMode', 'auto');
    saveas(gcf, fullFileName, 'svg')
end
