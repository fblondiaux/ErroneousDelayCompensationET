% Script to investigate the effect of inertia on the oscillations.
% Figure 5 of the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024

%% Simulation parameters
simParams;
nbForce = 3; % length(forces);

%% I = 0.15
I1 = 0.15;
%"Healthy Controls - HC"
delayError1 = 1; % Delay error in percentage
[freq, pert_x_HC, pert_xest_HC, pert_u_HC, PSD_HC] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I1, nbForce, x0, delayError1, delayError1);

%"Essential Tremor" - ET
delayError2 = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
[~, pert_x_ET, pert_xest_ET, pert_u_ET, PSD_ET] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I1, nbForce, x0, delayError2, delayError2);

%% I=0.10
I2 = 0.10;
%"Healthy Controls - HC"
[~, pert_x_HC10, pert_xest_HC10, pert_u_HC10, PSD_HC10] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I2, nbForce, x0, delayError1, delayError1);

%"Essential Tremor" - ET
[~, pert_x_ET10, pert_xest_ET10, pert_u_ET10, PSD_ET10] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I2, nbForce, x0, delayError2, delayError2);

%% I = 0.20
I3 = 0.20;
%"Healthy Controls - HC"
[~, pert_x_HC20, pert_xest_HC20, pert_u_HC20, PSD_HC20] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I3, nbForce, x0, delayError1, delayError1);

%"Essential Tremor" - ET
[~, pert_x_ET20, pert_xest_ET20, pert_u_ET20, PSD_ET20] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I3, nbForce, x0, delayError2, delayError2);


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
