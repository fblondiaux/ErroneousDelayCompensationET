% Calculate the discrete B matrix using a given A matrix, discretization step (dt), and B matrix.
% Parameters:
%   B: The B matrix.
%   A: The A matrix.
%   dt: The discretization step.
% Returns:
%   Bd: The discrete B matrix.
function Bd = discreteB(B, A, dt)
    % The routine uses a dt/100 discretization step
    % Correspond to eq. 9 of the paper
    step = dt / 100;
    delta = 0:step:dt;
    Bd = 0;

    for i = 1:size(delta, 2)
        Bd = Bd + step * expm(delta(i) * A);
    end

    Bd = Bd * B;
end
