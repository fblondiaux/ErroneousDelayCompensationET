% Script that investigates if the oscillations could be caused by an impaired
% state estimation where we underestimated either A or B. In these
% simulations, the delay is accurately estimated.
% Figure 5 in the publication.
%
% It is possible to run analyses for an under/over estimation of the Kalman
% gains using the same code. Change the value of the third element of the
% AB variable (last input of simulation). For example Kest = 0.2 or Kest =
% 2.5, use [1 1 0.2] or [1 1 2.5]
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

% A - Error only in A with Aest = 1/2A
scaleFactors = [0.5 1 1];
[~, pert_x_A, pert_xest_A, pert_u_A, PSD_A] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError, scaleFactors);

% B - Error only in B with Best = 1/2B
scaleFactors = [1 0.5 1];
[~, pert_x_B, pert_xest_B, pert_u_B, PSD_B] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError, scaleFactors);

%% Plot with the arm angle, angular velocity, control and PSD

tv = -10 * dt:dt:(timeStab - dt - 10 * dt); % Perturbation start after 10 timesteps
tv = tv * 1000; %From s to ms

constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332 * 86.11, 11/216 * 64.43);
ax = subplot(1, 11, 1:5, 'Units', 'centimeters');
ax.Position = [14.8, 66.5, 33.5, 23.86] / 10;
hold on;

for f = forces
    plot(tv, squeeze(mean(pert_x_HC(f, :, 1, :), 2)) * 180 / pi, 'Color', color_c, 'LineWidth', thickLine);
    plot(tv, squeeze(mean(pert_x_A(f, :, 1, :), 2)) * 180 / pi, 'Color', color_var1, 'LineWidth', thickLine); %,
    plot(tv, squeeze(mean(pert_x_B(f, :, 1, :), 2)) * 180 / pi, 'Color', color_var2, 'LineWidth', thickLine); %
end

xline(0);
xlabel('Time (ms)');
ylabel('Angle (deg)');
title('Elbow joint angle');
xlim([tv(1) tv(end)])

%For different error of the Kalman gains
%leg = legend(["No Error","K=0.2*K","K=2.5*K",],'FontSize',4, 'Location','best');

leg = legend(["No Error", "A=0.5*A", "B=0.5*B", ], 'FontSize', 4, 'Location', 'best');
leg.ItemTokenSize = [30/3, 18/3];
%figForInkscapeSave(F,append(figurePath,'runSimAngle'))

%% Power spectral density - Normalized
ax = subplot(1, 11, 6:7, 'Units', 'centimeters');
ax.Position = [56.3, 66.5, 15.2, 23.86] / 10; % define your position
hold on;

[M, I] = max(squeeze(mean(mean(PSD_A(:, :, :), 2), 1)));
disp('The peak of the frequency is at (Hz):')
disp(freq(I))

plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_A, 2), 1)) / M, 'Color', color_var1, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_B, 2), 1)) / M, 'Color', color_var2, 'LineWidth', thickLine)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('PSD');
hold on;
xlim([0 15]);

% Change the scale of the axis for values above 1. - For different errors of
% K
% yticks = 0:5;
% yticksLabel = [0 1 5 10 15 20]; % Final tick values
% YVal = squeeze(mean(mean(PSD_B,2),1))/M;
% YVal(YVal>1) = YVal(YVal>1)/5 +1;
%
% plot(freq,YVal,'Color',color_var2,'LineWidth',thickLine);
% ylim([0 max(yticks)]);
% set(gca,'YTick',yticks,'YTickLabel',yticksLabel); % Update the ticks

%% Control Signal R2
forces = 1:3;
%LLR - R2
ax = subplot(1, 11, 8:9, 'Units', 'centimeters');
ax.Position = [79.5, 66.5, 13, 23.86] / 10; % define your position
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_A(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_var1, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var1, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_B(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_var2, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var2, 'MarkerSize', 5);
axis([0.3 3.5 -5 25])
xlabel('Perturbation (Nm)')
ylabel('Control (a.u.)')
title('R2')

%LLR - R3
ax = subplot(1, 11, 10:11, 'Units', 'centimeters');
ax.Position = [100.5, 66.5, 13, 23.86] / 10; % define your position
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_A(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_var1, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var1, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_B(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_var2, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var2, 'MarkerSize', 5);
axis([0.5 3.5 -5 25])
title('R3')

sgtitle('Impact of an error in the internal prediction')
figForInkscapeSave(F, append(figurePath, 'allSim_AB'))
