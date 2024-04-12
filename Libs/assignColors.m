% Function to assign colors based on a given factor.
% Parameters:
%   the_force: The force value.
% Returns:
%   c: The assigned color for c.
%   p: The assigned color for p.
function [c, p] = assignColors(the_force)
    constantsPlots;
    if the_force == 1
        c = color_c_light;
        p = color_p_light;
    elseif the_force == 2
        c = color_c_medium;
        p = color_p_medium;
    else
        c = color_c;
        p = color_p;
    end
end