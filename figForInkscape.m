%Function used for making the figures, defines some parameters to the
%figure, dimension, units, ...

function [F] = figForInkscape(FigWidth,FigHeight)

F=figure('PaperPositionMode','manual');
AX = gca;
F.Resize = 'off';
F.Units = 'centimeters';
F.PaperUnits = 'centimeters';
AX.Units = 'centimeters';
AX.Clipping = 'on';
AX.PositionConstraint = 'innerposition';
Buffer = 5;
AX.InnerPosition = [Buffer Buffer FigWidth FigHeight]; % converting from pixels to points
F.OuterPosition = [0 0 FigWidth+2*Buffer FigHeight+2*Buffer]; % converting from pixels to points
F.PaperPosition = [0 0 FigWidth+2*Buffer FigHeight+2*Buffer];% converting from pixels to points

end