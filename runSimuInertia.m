% Script to investigate the effect of inertia on the oscillations.
% Figure 5 of the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024

%% Simulation parameters
nbSim = 20;
timeStab = 0.8; %Stabilization time
dt = 0.005; %[s]
delta = 0.055; %Delay of the feedback [s]
I = 0.15;
forces = [1 2 3]; % [Nm]
nbForce = 3; % length(forces);
x0 = [0; 0; 0; 0; 0; 0; 0; 0]; %Initial state vector, state variables include
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque

%% I = 0.15

%"Healthy Controls - HC"
pert_x_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC = zeros(nbForce, nbSim, 65);

for f = 1:3

    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est, xy, L_HC] = simulation(timeStab, delayError, delayError, f, x0, delta, I);
        pert_x_HC(f, i, :, :) = squeeze(x);
        pert_xest_HC(f, i, :, :) = squeeze(x_est);
        pert_u_HC(f, i, :, :) = squeeze(u);
        [freq, PSD_HC(f, i, :)] = getPSD(diff(pert_x_HC(f, i, 2, 40:end)) / dt, dt);
    end

end

%"Essential Tremor" - ET
pert_x_ET = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET = zeros(nbForce, nbSim, 65);
%With a error
for f = 1:3

    for i = 1:nbSim
        delayError = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
        [x, u, x_est, xy, L_ET] = simulation(timeStab, delayError, delayError, f, x0, delta, I);
        pert_x_ET(f, i, :, :) = squeeze(x);
        pert_xest_ET(f, i, :, :) = squeeze(x_est);
        pert_u_ET(f, i, :, :) = squeeze(u);
        [freq, PSD_ET(f, i, :)] = getPSD(diff(pert_x_ET(f, i, 2, 40:end)) / dt, dt);
    end

end

%% I=0.10
%"Healthy Controls - HC"
pert_x_HC10 = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC10 = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC10 = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC10 = zeros(nbForce, nbSim, 65);

for f = 1:3

    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est, xy, L_HC] = simulation(timeStab, delayError, delayError, f, x0, delta, .10);
        pert_x_HC10(f, i, :, :) = squeeze(x);
        pert_xest_HC10(f, i, :, :) = squeeze(x_est);
        pert_u_HC10(f, i, :, :) = squeeze(u);
        [freq, PSD_HC10(f, i, :)] = getPSD(diff(pert_x_HC10(f, i, 2, 40:end)) / dt, dt);
    end

end

%"Essential Tremor" - ET
pert_x_ET10 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET10 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET10 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET10 = zeros(nbForce, nbSim, 65);
%With a error
for f = 1:3

    for i = 1:nbSim
        delayError = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
        [x, u, x_est, xy, L_ET10] = simulation(timeStab, delayError, delayError, f, x0, delta, 0.10);
        pert_x_ET10(f, i, :, :) = squeeze(x);
        pert_xest_ET10(f, i, :, :) = squeeze(x_est);
        pert_u_ET10(f, i, :, :) = squeeze(u);
        [freq, PSD_ET10(f, i, :)] = getPSD(diff(pert_x_ET10(f, i, 2, 40:end)) / dt, dt);
    end

end

%% I = 0.20
%"Healthy Controls - HC"
pert_x_HC20 = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC20 = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC20 = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC20 = zeros(nbForce, nbSim, 65);

for f = 1:3

    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est, xy, L_HC20] = simulation(timeStab, delayError, delayError, f, x0, delta, .2);
        pert_x_HC20(f, i, :, :) = squeeze(x);
        pert_xest_HC20(f, i, :, :) = squeeze(x_est);
        pert_u_HC20(f, i, :, :) = squeeze(u);
        [freq, PSD_HC20(f, i, :)] = getPSD(diff(pert_x_HC20(f, i, 2, 40:end)) / dt, dt);
    end

end

%"Essential Tremor" - ET
pert_x_ET20 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET20 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET20 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET20 = zeros(nbForce, nbSim, 65);
%With a error
for f = 1:3

    for i = 1:nbSim
        delayError = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
        [x, u, x_est, xy, L_ET20] = simulation(timeStab, delayError, delayError, f, x0, delta, 0.2);
        pert_x_ET20(f, i, :, :) = squeeze(x);
        pert_xest_ET20(f, i, :, :) = squeeze(x_est);
        pert_u_ET20(f, i, :, :) = squeeze(u);
        [freq, PSD_ET20(f, i, :)] = getPSD(diff(pert_x_ET20(f, i, 2, 40:end)) / dt, dt);
    end

end

%% Plot the PSD
constantsPlots;
F = figForInkscape(19/332 * 86.11, 10/216 * 64.43);

ax = subplot(2, 10, 1:20, 'Units', 'centimeters');
ax.Position = [14.8, 26, 25.2, 34.4] / 10; % define your position
hold on;

% [~,I]= max(squeeze(mean(mean(PSD_ET10,2),1)));
% %disp('The peak of the frequency is at (Hz):');disp(freq(I))
%
% [~,I]= max(squeeze(mean(mean(PSD_ET20,2),1)));
% %disp('The peak of the frequency is at (Hz):');disp(freq(I))

[M, I] = max(squeeze(mean(mean(PSD_ET, 2), 1)));
%disp('The peak of the frequency is at (Hz):')
%disp(freq(I))

%Dummy print in black for the legend
plot15 = plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine);
plot10 = plot(freq, squeeze(mean(mean(PSD_HC10, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine, 'Linestyle', "-.");
plot20 = plot(freq, squeeze(mean(mean(PSD_HC20, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine, 'Linestyle', "--");

%Print in the actual colors
plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine);
plot(freq, squeeze(mean(mean(PSD_ET, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine)

plot(freq, squeeze(mean(mean(PSD_HC10, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine, 'Linestyle', "-.");
plot(freq, squeeze(mean(mean(PSD_ET10, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine, 'Linestyle', "-.")

plot(freq, squeeze(mean(mean(PSD_HC20, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine, 'Linestyle', "--");
plot(freq, squeeze(mean(mean(PSD_ET20, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine, 'Linestyle', "--")

xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Variation of Inertia');
hold on;
xlim([0 15]);
leg = legend([plot10, plot15, plot20], {'I=0.1', 'I=0.15', 'I=2'});
leg.ItemTokenSize = [30/2, 18/2];
figForInkscapeSave(F, append(figurePath, 'allSim_Inertia'))
