% Script to investigate the effect of delay duration on the oscillations.
% Figure 5 of the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024
%

%% Simulation parameters
simParams;
nbForce = 3; % length(forces);

%% Delay 55 ms
delta1 = 0.055; %Delay of the feedback [s]
%"Healthy Controls - HC"
delayError1 = 1; % Delay error in percentage
[freq, pert_x_HC, pert_xest_HC, pert_u_HC, PSD_HC] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta1, I, nbForce, x0, delayError1, delayError1);


%"Essential Tremor" - ET
delayError2 = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
[~, pert_x_ET, pert_xest_ET, pert_u_ET, PSD_ET] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta1, I, nbForce, x0, delayError2, delayError2);

%% Delay 65 ms
delta2 = 0.065; %Delay of the feedback [s]
%"Healthy Controls - HC"
[~, pert_x_H65, pert_xest_HC65, pert_u_HC65, PSD_HC65] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta2, I, nbForce, x0, delayError1, delayError1);

%"Essential Tremor" - ET
[~, pert_x_ET65, pert_xest_ET65, pert_u_ET65, PSD_ET65] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta2, I, nbForce, x0, delayError2, delayError2);

%% Delay 45 ms
delta3 = 0.045; %Delay of the feedback [s]
%"Healthy Controls - HC"
[~, pert_x_HC45, pert_xest_HC45, pert_u_HC45, PSD_HC45] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta3, I, nbForce, x0, delayError1, delayError1);

%"Essential Tremor" - ET
[~, pert_x_ET45, pert_xest_ET45, pert_u_ET45, PSD_ET45] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta3, I, nbForce, x0, delayError2, delayError2);

%% Plot with the arm angle, angular velocity, control and PSD
%% PSD
constantsPlots;
F = figForInkscape(19/332 * 86.11, 10/216 * 64.43);

ax = subplot(2, 10, 1:20, 'Units', 'centimeters');
ax.Position = [14.8, 26, 25.2, 34.4] / 10; % define your position
hold on;

% [M,I]= max(squeeze(mean(mean(PSD_ET65,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
%
% [M,I]= max(squeeze(mean(mean(PSD_ET45,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))

[M, I] = max(squeeze(mean(mean(PSD_ET, 2), 1)));
% disp('The peak of the frequency is at (Hz):')
% disp(freq(I))

% Dummy plots to have the legend in black
plot45 = plot(freq, squeeze(mean(mean(PSD_HC45, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine, 'Linestyle', "--");
plot55 = plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine);
plot65 = plot(freq, squeeze(mean(mean(PSD_HC65, 2), 1)) / M, 'Color', 'k', 'LineWidth', thickLine, 'Linestyle', "-.");

%Actual plots
plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine);
plot(freq, squeeze(mean(mean(PSD_ET, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine)

plot(freq, squeeze(mean(mean(PSD_HC65, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine, 'Linestyle', "-.");
plot(freq, squeeze(mean(mean(PSD_ET65, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine, 'Linestyle', "-.")

plot(freq, squeeze(mean(mean(PSD_HC45, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine, 'Linestyle', "--");
plot(freq, squeeze(mean(mean(PSD_ET45, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine, 'Linestyle', "--")

xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Variation of delay duration');
hold on;
xlim([0 15]);
leg = legend([plot45, plot55, plot65], {'d=45ms', 'd=55ms', 'd=65ms'});
leg.ItemTokenSize = [30/2, 18/2];
figForInkscapeSave(F, append(figurePath, 'allSim_DelayDuration'))
