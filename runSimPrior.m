% Script that investigates if the oscillations could be caused by an impaired
% state estimation where we underestimated or overestimated the prior. In these
% simulations, the delay is accurately estimated.
% Figure 5 in the publication.
%
% It is possible to run analyses for an under/over estimation of the Kalman
% gains using the same code. Change the value of the second element of the
% errorEstimation variable (last input of simulation). For example Kest = 0.2 or Kest =
% 2.5, use [1 0.2] or [1 2.5]
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

% Underestimated Prior
scaleFactors = [0.7 1];
[~, pert_x_underEst, pert_xest_underEst, pert_u_underEst, PSD_underEst] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError, scaleFactors);

% Overestimated prior
scaleFactors = [1.1 1];
[~, pert_x_overEst, pert_xest_overEst, pert_u_overEst, PSD_overEst] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError, scaleFactors);

%% Plot with the arm angle, angular velocity, control and PSD

tv = -10 * dt:dt:(timeStab - dt - 10 * dt); % Perturbation start after 10 timesteps
tv = tv * 1000; %From s to ms

constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332 * 86.11, 11/216 * 64.43);

% For angle plot
[ax1] = plotData(forces, tv, pert_x_HC, pert_x_underEst, 1, 'Angle (deg)', 'Elbow joint angle', ...
    thickLine, [14.8, 66.5, 33.5, 23.86] / 10, {2, 11, 1:5}, ...
    color_c, color_var1, color_var2, pert_x_overEst);

%For different error of the Kalman gains
%leg = legend(["No Error","K=0.2*K","K=2.5*K",],'FontSize',4, 'Location','best');

leg = legend(["No Error","prior=0.7*prior","prior=1.1*prior",],'FontSize',4, 'Location','best');
leg.ItemTokenSize = [30/3, 18/3];
%figForInkscapeSave(F,append(figurePath,'runSimAngle'))

%% Power spectral density - Normalized
ax = subplot(1, 11, 6:7, 'Units', 'centimeters');
ax.Position = [56.3, 66.5, 15.2, 23.86] / 10; % define your position
hold on;

[M, I] = max(squeeze(mean(mean(PSD_underEst(:, :, :), 2), 1)));
disp('The peak of the frequency is at (Hz):')
disp(freq(I))

plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_underEst, 2), 1)) / M, 'Color', color_var1, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_overEst, 2), 1)) / M, 'Color', color_var2, 'LineWidth', thickLine)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('PSD');
hold on;
xlim([0 15]);

% Change the scale of the axis for values above 1. - For different errors of
% K
% yticks = 0:5;
% yticksLabel = [0 1 5 10 15 20]; % Final tick values
% YVal = squeeze(mean(mean(PSD_overEst,2),1))/M;
% YVal(YVal>1) = YVal(YVal>1)/5 +1;
%
% plot(freq,YVal,'Color',color_var2,'LineWidth',thickLine);
% ylim([0 max(yticks)]);
% set(gca,'YTick',yticks,'YTickLabel',yticksLabel); % Update the ticks

%% Control Signal R2
forces = 1:3;
%LLR - R2
% For R2 plot
plotErrorBar(forces, pert_u_HC, color_c, pert_u_underEst, color_var1, thickLine,...
 [79.5, 66.5, 13, 23.86] / 10, {1, 11, 8:9}, 'Control (a.u.)', 'R2', 19:25, [0.3 3.5 -5 25], pert_u_overEst, color_var2);

% For R3 plot
plotErrorBar(forces, pert_u_HC, color_c, pert_u_underEst, color_var1, thickLine,...
 [100.5, 66.5, 13, 23.86] / 10, {1, 11, 10:11}, 'Control (a.u.)', 'R3', 25:31, [0.5 3.5 -5 25], pert_u_overEst, color_var2);

%% Title

sgtitle('Impact of an error in the internal prediction')

%Save the figure
savefigure(F, figurePath, 'allSim_prior');