% Integration of the state and motor commands over the delay interval
% Corresponds to eq. 11 of the paper.
%
% Parameters:
%   - y = sensory feedback
%   - dt = discretization step
%   - delay  = duration of the delay
%   - u = motor commands
%   - A, B = System dynamics
%   - pErr1 = proportional error in the delay parameter - integration of the state
%   - pErr2 = proportional error in the delay parameter - integration of the motor commands
%
% Returns:
%   - xy = integrated state
%   - V = integrated variance
%   - M = exponential matrix
function [xy, V, M] = integrate(y, dt, delay, u, A, B, pErr1, pErr2)

    %If want to separate the contribution of a delay error on the
    %integration of the state from the integration of the motor command
    delay1 = pErr1 * delay; % Delay used for M*y
    delay2 = pErr2 * delay; % Delay used to integrate motor commands

    %% ALT: Uncomment one of the following line to add noise to the delay estimation.
    %noise = [0.005 0.010 0 -0.005 -0.010];
    %noise = [0.005 0.010 0.015 0 -0.005 -0.010 -0.015];
    noise = 0;
    r = randi([1 length(noise)]);
    delay1 = delay1 + noise(r);
    delay2 = delay2 + noise(r);

    V = 0;
    xy = 0;
    count = 0;
    M = expm(delay1 * A);

    delay_d = round(delay2 / dt); % [nb Steps]

    if length(u) < delay_d + 1
        % Deal with cases where there isn't enough motor commands (beginning of the movement)
        temp_u = ones(1, delay_d + 1) * u(1);
        temp_u(end - length(u) + 1:end) = u;
        u = temp_u;
    end

    for s = 0:dt:delay2

        count = count + 1;
        count = max(1, count);
        count = min(count, size(u, 2));
        integrant = M * expm(-s * A) * B * u(count);
        xy = xy + dt * integrant;
        V = V + dt * (integrant * integrant');

    end

    xy = xy + M * y;
