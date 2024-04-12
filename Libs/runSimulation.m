% Run simulation with specified parameters and return the frequency, perturbations, and power spectral density.
%
% Parameters:
%   - nbSim: number of simulations
%   - nbState: number of states
%   - nbControl: number of controls
%   - timeStab: stabilization time
%   - dt: time step between samples in seconds
%   - delta: Delay of the feedback in seconds
%   - I: Inertia value in Kgm²
%   - nbForce: number of forces
%   - x0: initial state vector
%   - delayError1: delay error value 1
%   - delayError2: delay error value 2
%   - errorEstimation = scaling factors to create mismatch in the estimation of the prior 
%     and kalman gains, [1 1] = no errors on the estimation.
% Returns:
%   - freq: frequency of the simulation
%   - pert_x: state vector x
%   - pert_xest: estimated state vector x_est
%   - pert_u: control vector u
%   - PSD: power spectral density of the acceleration
function [freq, pert_x, pert_xest, pert_u, PSD, xy, L] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError1, delayError2, scaleFactors)
    pert_x = zeros(nbForce, nbSim, nbState, round(timeStab / dt));
    pert_xest = zeros(nbForce, nbSim, nbState, round(timeStab / dt));
    pert_u = zeros(nbForce, nbSim, nbControl, round(timeStab / dt));
    PSD = zeros(nbForce, nbSim, nbState * nbState + nbControl);

    % Check if errorEstimation is defined
    if nargin == 11
        scaleFactors = [1 1];
    end

    for f = 1:nbForce

        for i = 1:nbSim
            [x, u, x_est, xy, L] = simulation(timeStab, delayError1, delayError2, f, x0, delta, I, scaleFactors);
            pert_x(f, i, :, :) = squeeze(x);
            pert_xest(f, i, :, :) = squeeze(x_est);
            pert_u(f, i, :, :) = squeeze(u);
            [freq, PSD(f, i, :)] = getPSD(diff(pert_x(f, i, 2, 40:end)) / dt, dt); % PSD of the angular velocity
        end

    end

end
