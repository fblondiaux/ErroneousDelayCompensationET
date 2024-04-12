% Main script for generating the simulations
%
% Author: Flo Blondiaux
% Date: Jan 2024
%
% This code performs simulations for a postural perturbation task, for several
% perturbation magnitude with two conditions, accurate or under estimation
% of the delay in the state estimator.
% This code is used to generate the figure 4 of the paper.


% ALT COMMENTS
%Extreme changes in the slopes of the R3 response were obtained
% changing the following simulation parameters (Figure 4 - Dotted lines):
%Smaller slope
%I = 0.2;
%delta = 0.045;
%Bigger slope
%I=0.1;
%delta = 0.065;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation parameters
simParams;
nbForce = length(forces);

% "Healthy Controls - HC"
delayError = 1; % Delay error in percentage
[freq, pert_x_HC, pert_xest_HC, pert_u_HC, PSD_HC] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError);

% "Essential Tremor" - ET
delayError = .7; % Delay is underestimated : Delay used in state estimation is 70% of the actual delay
[~, pert_x_ET, pert_xest_ET, pert_u_ET, PSD_ET] = runSimulation(nbSim, nbState, nbControl, timeStab, dt, delta, I, nbForce, x0, delayError, delayError);

%% Plot with the arm angle, angular velocity, control and PSD

tv = -10 * dt:dt:(timeStab - dt - 10 * dt); % Perturbation start after 10 timesteps
tv = tv * 1000; %From s to ms
constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332 * 86.11, 10/216 * 64.43);
ax = subplot(2, 11, 1:5, 'Units', 'centimeters');
ax.Position = [14.8, 66.5, 33.5, 23.86] / 10; % define your position
hold on;

for f = forces

    if f == 1
        col_c = color_c_light;
        col_p = color_p_light;
    elseif f == 2
        col_c = color_c_medium;
        col_p = color_p_medium;
    else
        col_c = color_c;
        col_p = color_p;
    end

    plot(tv, squeeze(mean(pert_x_HC(f, :, 1, :), 2)) * 180 / pi, 'Color', col_c, 'LineWidth', thickLine); %From rad to °
    plot(tv, squeeze(mean(pert_x_ET(f, :, 1, :), 2)) * 180 / pi, 'Color', col_p, 'LineWidth', thickLine);
    xline(0);
    xlabel('Time (ms)');
    ylabel('Angle (deg)');
    title('Elbow joint angle');
    xlim([tv(1) tv(end)])
end

%% Plot of the velocity
ax = subplot(2, 11, 6:10, 'Units', 'centimeters');
ax.Position = [61, 66.5, 33.5, 23.86] / 10;
hold on

for f = forces

    if f == 1
        col_c = color_c_light;
        col_p = color_p_light;
    elseif f == 2
        col_c = color_c_medium;
        col_p = color_p_medium;
    else
        col_c = color_c;
        col_p = color_p;
    end

    plot(tv, squeeze(mean(pert_x_HC(f, :, 2, :), 2)) * 180 / pi, 'Color', col_c, 'LineWidth', thickLine);
    plot(tv, squeeze(mean(pert_x_ET(f, :, 2, :), 2)) * 180 / pi, 'Color', col_p, 'LineWidth', thickLine);
    xline(0);
    xlabel('Time (ms)');
    ylabel('Velocity (deg/s)');
    title('Elbow joint velocity');
    xlim([tv(1) tv(end)])
end

%% Power spectral density - Normalized
ax = subplot(2, 11, 12:13, 'Units', 'centimeters'); % define your position
ax.Position = [14.8, 26, 15.2, 24.4] / 10;
hold on;

[M, I] = max(squeeze(mean(mean(PSD_ET(:, :, :), 2), 1))); %Get the peak value for normalization
% disp('The peak of the frequency is at (Hz):')
% disp(freq(I))

plot(freq, squeeze(mean(mean(PSD_HC, 2), 1)) / M, 'Color', color_c, 'LineWidth', thickLine)
plot(freq, squeeze(mean(mean(PSD_ET, 2), 1)) / M, 'Color', color_p, 'LineWidth', thickLine)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Simulations');
hold on;
xlim([0 15]); %Limit to 15 Hz

%% Control signal
ax = subplot(2, 11, 14:17, 'Units', 'centimeters');
ax.Position = [38, 26, 24.4, 24.4] / 10;

for f = 1:3

    if f == 1
        col_c = color_c_light;
        col_p = color_p_light;
    elseif f == 2
        col_c = color_c_medium;
        col_p = color_p_medium;
    else
        col_c = color_c;
        col_p = color_p;
    end

    plot(tv(1:40), squeeze(mean(pert_u_HC(f, :, 1, 1:40), 2)), 'Color', col_c, 'LineWidth', thickLine);
    hold on;
    plot(tv(1:40), squeeze(mean(pert_u_ET(f, :, 1, 1:40), 2)), 'Color', col_p, 'LineWidth', thickLine);
end

xlabel('Time (ms)');
ylabel('Control (a.u.)');
title('Control Signal');
xline(0)
xline(25)
xline(45)
xline(75)
xline(105)
ylim([-0.5 20]);

%% LLR - R2
ax = subplot(2, 11, 18:19, 'Units', 'centimeters');
ax.Position = [70.4, 26, 13, 24.4] / 10;
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_ET(:, :, 1, 19:25), 4), 2))', [0 0 0], '-o', 'Color', color_p, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_p, 'MarkerSize', 5);
axis([0.3 3.5 -5 25])
xlabel('Perturbation (Nm)')
ylabel('Control (a.u.)')
title('R2')

%% LLR - R3
ax = subplot(2, 11, 20:21, 'Units', 'centimeters');
ax.Position = [88, 26, 13, 24.4] / 10;
hold on;
errorbar(forces, squeeze(mean(mean(pert_u_HC(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_c, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_c, 'MarkerSize', 5);
errorbar(forces, squeeze(mean(mean(pert_u_ET(:, :, 1, 25:31), 4), 2))', [0 0 0], '-o', 'Color', color_p, 'LineWidth', thickLine, 'CapSize', 0, 'MarkerFaceColor', color_p, 'MarkerSize', 5);
axis([0.5 3.5 -5 25])
title('R3')

%% Estimated state
%Uncomment to add small pannels on the side with the estimted state over R2
%and R3 Windows

% ax = subplot(2,11,11,'Units','centimeters');
% ax.Position =[61+60,60+20,33.5*0.35,23.86*0.35]/10;
% hold on;
% plot(tv(20:32),squeeze(mean(pert_xest_ET(f, :, 1, 19:31), 2))*180/pi,'Color',color_p,'LineWidth',thickLine); %,
% plot(tv(20:32),squeeze(mean(pert_xest_HC(f, :, 1, 19:31),2))*180/pi,'Color',color_c,'LineWidth',thickLine);
%
% ax = subplot(2,11,22,'Units','centimeters'); % define your position
% ax.Position =[61+60,20,33.5*0.35,23.86*0.35]/10;
% hold on;
% plot(tv(20:32),squeeze(mean(pert_xest_ET(f, :, 2, 19:31), 2))*180/pi,'Color',color_p,'LineWidth',thickLine); %,
% plot(tv(20:32),squeeze(mean(pert_xest_HC(f, :, 2, 19:31),2))*180/pi,'Color',color_c,'LineWidth',thickLine);
% ylim([-100 0])
% xlim([45 105])
%
%check if we are on Mac or Windows
figname = 'allSim';

if ismac
    figname = 'allSim.fig';
end

figForInkscapeSave(F, append(figurePath, figname))
