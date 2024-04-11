% Script that investigates the impact of the origin of the delay error
% Either solely in the integration of the state over the delay interval
% Either solely in the integration of the motor commands over the delay
% interval
% Figure 5 in the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024

%% Simulation parameters
nbSim = 20;
timeStab = 0.8; %Stabilization time
dt = 0.005; %[s]
delta = 0.055; %Delay of the feedback [s]
I = 0.15;
forces = [2]; % [Nm]
nbForce = 3; % length(forces);
x0 = [0; 0; 0; 0; 0; 0; 0; 0]; %Initial state vector, state variables include
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque

%"Healthy Controls - HC"
pert_x_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC = zeros(nbForce, nbSim, 1, round((timeStab) / dt)); %1 = nb of control variables
PSD_HC = zeros(nbForce, nbSim, 65);

for f = 1:3

    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est] = simulation(timeStab, delayError, delayError, f, x0, delta, I);
        pert_x_HC(f, i, :, :) = squeeze(x);
        pert_xest_HC(f, i, :, :) = squeeze(x_est);
        pert_u_HC(f, i, :, :) = squeeze(u);
        [freq, PSD_HC(f, i, :)] = getPSD(diff(pert_x_HC(f, i, 2, 40:end)) / dt, dt); % PSD computed on the stabilization time
    end

end

% Y - Error only in M*y
pert_x_Y = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_xest_Y = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_Y = zeros(nbForce, nbSim, 1, round((timeStab) / dt));
PSD_Y = zeros(nbForce, nbSim, 65);
%With a error
for f = 1:3

    for i = 1:nbSim
        delayError = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
        [x, u, x_est] = simulation(timeStab, 0.7, 1, f, x0, delta, I);
        pert_x_Y(f, i, :, :) = squeeze(x);
        pert_xest_Y(f, i, :, :) = squeeze(x_est);
        pert_u_Y(f, i, :, :) = squeeze(u);
        [freq, PSD_Y(f, i, :)] = getPSD(diff(pert_x_Y(f, i, 2, 40:end)) / dt, dt);
    end

end

% Sig- error only in the Sum M*B*u
pert_x_Sig = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_xest_Sig = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_Sig = zeros(nbForce, nbSim, 1, round((timeStab) / dt));
PSD_Sig = zeros(nbForce, nbSim, 65);
%With a error
for f = 1:3

    for i = 1:nbSim
        delayError = .7; % Delay is underestimated : Delay used is 70 % of the actual delay
        [x, u, x_est] = simulation(timeStab, 1, 0.7, f, x0, delta, I);
        pert_x_Sig(f, i, :, :) = squeeze(x);
        pert_xest_Sig(f, i, :, :) = squeeze(x_est);
        pert_u_Sig(f, i, :, :) = squeeze(u);
        [freq, PSD_Sig(f, i, :)] = getPSD(diff(pert_x_Sig(f, i, 2, 40:end)) / dt, dt);
    end

end

%% Plot with the arm angle, angular velocity, control and PSD

tv = -10 * dt:dt:(timeStab - dt - 10 * dt); % Perturbation start after 10 timesteps
tv = tv * 1000; %From s to ms

constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332 * 86.11, 11/216 * 64.43);
ax = subplot(1, 11, 1:5, 'Units', 'centimeters');
ax.Position = [14.8, 66.5, 33.5, 23.86] / 10; % define your position
hold on;

for f = forces
    plot(tv, squeeze(mean(pert_x_HC(f, :, 1, :), 2)) * 180 / pi, 'Color', color_c, 'LineWidth', thickLine);
    plot(tv, squeeze(mean(pert_x_Y(f, :, 1, :), 2)) * 180 / pi, 'Color', color_var1, 'LineWidth', thickLine); %,
    plot(tv, squeeze(mean(pert_x_Sig(f, :, 1, :), 2)) * 180 / pi, 'Color', color_var2, 'LineWidth', thickLine); %

end

xline(0);
xlabel('Time (ms)');
ylabel('Angle (deg)');
title('Elbow joint angle');
xlim([tv(1) tv(end)])
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
%LLR - R2
ax = subplot(1, 11, 8:9, 'Units', 'centimeters');
ax.Position = [79.5, 66.5, 13, 23.86] / 10; % define your position
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_Y(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_var1, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var1, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_Sig(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_var2, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var2, 'MarkerSize', 5);
axis([0.3 3.5 -5 25])
xlabel('Perturbation (Nm)')
ylabel('Control (a.u.)')
title('R2')

%LLR - R3
ax = subplot(1, 11, 10:11, 'Units', 'centimeters');
ax.Position = [100.5, 66.5, 13, 23.86] / 10; % define your position
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_Y(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_var1, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var1, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_Sig(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_var2, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_var2, 'MarkerSize', 5);
axis([0.5 3.5 -5 25])
title('R3')

sgtitle('Impact of the origin of the delay compensation error')
figForInkscapeSave(F, append(figurePath, 'allSim_2Pathways'))
