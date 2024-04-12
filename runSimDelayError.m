% Script to investigate the effect of delay error on the oscillations.
% Figure 5 of the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024


%% Simulation parameters
simParams;
nbForce = 3; % length(forces);

%"Healthy Controls - HC"
delayError1 = 1; % Delay error in percentage
[freq, pert_x_HC, pert_xest_HC, pert_u_HC, PSD_HC] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError1, delayError1);


%"Essential Tremor" - ET
%% Delay error = 90%
delayError2 = .9; % Delay is underestimated : Delay used is 90 % of the actual delay
[~, pert_x_ET9, pert_xest_ET9, pert_u_ET9, PSD_ET9] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError2, delayError2);

%% Delay error = 80%
delayError3 = .8; % Delay is underestimated : Delay used is 80 % of the actual delay
[~, pert_x_ET8, pert_xest_ET8, pert_u_ET8, PSD_ET8] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError3, delayError3);

%% Delay error = 70%
delayError4 = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
[~, pert_x_ET7, pert_xest_ET7, pert_u_ET7, PSD_ET7] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError4, delayError4);

%% Delay error = 60%
delayError5 = .6; % Delay is underestimated : Delay used is 60 % of the actual delay
[~, pert_x_ET6, pert_xest_ET6, pert_u_ET6, PSD_ET6] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError5, delayError5);

%% Plot with the arm angle, angular velocity, control and PSD

%% Power spectral density - Normalized
constantsPlots;
F = figForInkscape(19/332 * 86.11, 10/216 * 64.43);
ax = subplot(2, 10, 1:20, 'Units', 'centimeters');
ax.Position = [14.8, 26, 25.2, 34.4] / 10; % define your position
hold on;
% [~,I]= max(squeeze(mean(mean(PSD_ET8,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
% [~,I]= max(squeeze(mean(mean(PSD_ET9,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
% [~,I]= max(squeeze(mean(mean(PSD_ET6,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))

[M, I] = max(squeeze(mean(mean(PSD_ET7, 2), 1)));
%disp('The peak of the frequency is at (Hz):')
%disp(freq(I))

plot(freq, squeeze(mean(mean(PSD_ET6, 2), 1)) / M, 'Color', color_var2, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_ET7, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_ET8, 2), 1)) / M, 'Color', color_var3, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_ET9, 2), 1)) / M, 'Color', color_var1, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine)
leg = legend(["^d = 0.6d", "^d = 0.7d", "^d = 0.8d", "^d = 0.9d", "^d = 1d"], 'FontSize', 4, 'Location', 'best');
leg.ItemTokenSize = [30/3, 18/3];
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Variation of delay errors');
hold on;
xlim([0 15]);

%Save the figure
savefigure(F, figurePath, 'allSim_delayErrors');