% Save the figure with the given name and path.
% Parameters:
% - F: the figure handle
% - figPath: the path to save the figure
% - figname: the name of the figure
function savefigure(F, figPath, figname)
    % check if we are on Mac or Windows
    if ismac
        figname = append(figname, '.fig');
    end

    figForInkscapeSave(F, append(figPath, figname))
end
