% Function used for making the figures, defines some parameters to the
% figure, dimension, units, ...
%
% Parameters:
%   - FigWidth: Width of the figure in centimeters.
%   - FigHeight: Height of the figure in centimeters.
%
% Returns:
%   - F: The created figure object.
function [F] = figForInkscape(FigWidth, FigHeight)
    
    % Create the figure and axes objects
    F = figure('PaperPositionMode', 'manual');
    AX = gca;
    F.Resize = 'off';
    F.Units = 'centimeters';
    F.PaperUnits = 'centimeters';
    AX.Units = 'centimeters';
    AX.Clipping = 'on';
    AX.PositionConstraint = 'innerposition';
    Buffer = 5;

    % converting from pixels to points
    %% 1 point = 1/72 inch and 1 inch = 2.54 cm
    %% This vector defines the size and position of the axes within the figure window
    AX.InnerPosition = [Buffer Buffer FigWidth FigHeight];
    F.OuterPosition = [0 0 FigWidth + 2 * Buffer FigHeight + 2 * Buffer];
    F.PaperPosition = [0 0 FigWidth + 2 * Buffer FigHeight + 2 * Buffer]; 

end
