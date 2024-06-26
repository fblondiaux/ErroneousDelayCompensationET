% Script that investigates the impact of the origin of the delay error
% Either solely in the integration of the state over the delay interval
% Either solely in the integration of the motor commands over the delay
% interval
% Figure 5 in the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024

%% Simulation parameters
simParams;
forces = [2]; % [Nm]
nbForce = 3; % length(forces);

%"Healthy Controls - HC"
delayError = 1; % Delay error in percentage
[freq, pert_x_HC, pert_xest_HC, pert_u_HC, PSD_HC] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError);

% Y - Error only in M*y
delayError1 = 0.7; % Delay is underestimated : Delay used is 70 % of the actual delay
delayError2 = 1; % Delay error in percentage
[~, pert_x_Y, pert_xest_Y, pert_u_Y, PSD_Y] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError1, delayError2);

% Sig- error only in the Sum M*B*u
delayError1 = 1; % Delay is underestimated : Delay used is 70 % of the actual delay
delayError2 = 0.7; % Delay error in percentage
[~, pert_x_Sig, pert_xest_Sig, pert_u_Sig, PSD_Sig] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError1, delayError2);

%% Plot with the arm angle, angular velocity, control and PSD

tv = -10 * dt:dt:(timeStab - dt - 10 * dt); % Perturbation start after 10 timesteps
tv = tv * 1000; %From s to ms

constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332 * 86.11, 11/216 * 64.43);

% For angle plot
[ax1] = plotData(forces, tv, pert_x_HC, pert_x_Y, 1, 'Angle (deg)', 'Elbow joint angle', ...
    thickLine, [14.8, 66.5, 33.5, 23.86] / 10, {2, 11, 1:5}, ...
     color_c, color_var1, color_var2, pert_x_Sig);

leg = legend(["No Error", "Extrap. y", "Extrap. u", ], 'FontSize', 4, 'Location', 'best');
leg.ItemTokenSize = [30/3, 18/3];

%% Power spectral density - Normalized
ax = subplot(1, 11, 6:7, 'Units', 'centimeters');
ax.Position = [56.3, 66.5, 15.2, 23.86] / 10; % define your position
hold on;

[M, I] = max(squeeze(mean(mean(PSD_Y(:, :, :), 2), 1)));
disp('The peak of the frequency is at (Hz):')
disp(freq(I))

plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_Y, 2), 1)) / M, 'Color', color_var1, 'LineWidth', thickLine)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('PSD');
hold on;
xlim([0 15]);

%Change the scale of the axis for values above 1.
yticks = 0:4;
yticksLabel = [0 1 10 20 30]; % Final tick values
YVal = squeeze(mean(mean(PSD_Sig, 2), 1)) / M;
YVal(YVal > 1) = YVal(YVal > 1) / 10 +1;

plot(freq, YVal, 'Color', color_var2, 'LineWidth', thickLine);
ylim([0 max(yticks)]);
set(gca, 'YTick', yticks, 'YTickLabel', yticksLabel); % Update the ticks

%% Control Signal R2
forces = 1:3;

% For R2 plot
plotErrorBar(forces, pert_u_HC, color_c, pert_u_Y, color_var1, thickLine,...
 [79.5, 66.5, 13, 23.86] / 10, {1, 11, 8:9}, 'Control (a.u.)', 'R2', 19:25, [0.2 3.5 -5 25], pert_u_Sig, color_var2);

% For R3 plot
plotErrorBar(forces, pert_u_HC, color_c, pert_u_Y, color_var1, thickLine,...
 [100.5, 66.5, 13, 23.86] / 10, {1, 11, 10:11}, 'Control (a.u.)', 'R3', 25:31, [0.5 3.5 -5 25], pert_u_Sig, color_var2);

%% Title
sgtitle('Impact of the origin of the delay compensation error')

%Save the figure
savefigure(F, figurePath, 'allSim_2Pathways');
