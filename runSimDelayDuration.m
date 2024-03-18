% Script to investigate the effect of delay duration on the oscillations. 
% Figure 5 of the publication.
%
% Author: Flo Blondiaux
% Date: Jan 2024
%

%% Simulation parameters
nbSim = 20;
timeStab = 0.8; %Stabilization time
dt = 0.005; %[s]
delta = 0.055; %Delay of the feedback [s]
I = 0.15; 
forces = [1 2 3]; % [Nm]
nbForce = 3;% length(forces);
x0 = [0;0;0;0;0;0;0;0]; %Initial state vector, state variables include 
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque


%% Delay 55 ms
%"Healthy Controls - HC"
pert_x_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC = zeros(nbForce,nbSim,65);

for f = 1:3
    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est] = simulation(timeStab, delayError, delayError, f, x0,delta,I);
        pert_x_HC(f, i, :, :) = squeeze(x);
        pert_xest_HC(f, i, :, :) = squeeze(x_est);
        pert_u_HC(f, i, :, :) = squeeze(u);
        [freq,PSD_HC(f,i,:)]= getPSD(diff(pert_x_HC(f, i, 2, 40:end))/dt,dt);
    end
end



%"Essential Tremor" - ET
pert_x_ET = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.7; % Delay is underestimated : Delay used is 70% of the actual delay
        [x, u, x_est] = simulation(timeStab, delayError,delayError, f, x0,delta,I);
        pert_x_ET(f, i, :, :) = squeeze(x);
        pert_xest_ET(f, i, :, :) = squeeze(x_est);
        pert_u_ET(f, i, :, :) = squeeze(u);
        [freq,PSD_ET(f,i,:)]= getPSD(diff(pert_x_ET(f, i, 2, 40:end))/dt,dt);
    end
end

%% Delay 65 ms
%"Healthy Controls - HC"
pert_x_H65 = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC65 = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC65 = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC65 = zeros(nbForce,nbSim,65);

for f = 1:3
    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est] = simulation(timeStab, delayError, delayError,f, x0,0.065,I);
        pert_x_H65(f, i, :, :) = squeeze(x);
        pert_xest_HC65(f, i, :, :) = squeeze(x_est);
        pert_u_HC65(f, i, :, :) = squeeze(u);
        [freq,PSD_HC65(f,i,:)]= getPSD(diff(pert_x_H65(f, i, 2, 40:end))/dt,dt);
    end
end


%"Essential Tremor" - ET
pert_x_ET65 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET65 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET65 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET65 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.7; % Delay is underestimated : Delay used is 70% of the actual delay
        [x, u, x_est] = simulation(timeStab, delayError,delayError, f, x0,0.065,I);
        pert_x_ET65(f, i, :, :) = squeeze(x);
        pert_xest_ET65(f, i, :, :) = squeeze(x_est);
        pert_u_ET65(f, i, :, :) = squeeze(u);
        [freq,PSD_ET65(f,i,:)]= getPSD(diff(pert_x_ET65(f, i, 2, 40:end))/dt,dt);
    end
end

%% Delay 45 ms
%"Healthy Controls - HC"
pert_x_HC45 = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC45 = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC45 = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC45 = zeros(nbForce,nbSim,65);

for f = 1:3
    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est] = simulation(timeStab, delayError, delayError,f, x0,0.045,I);
        pert_x_HC45(f, i, :, :) = squeeze(x);
        pert_xest_HC45(f, i, :, :) = squeeze(x_est);
        pert_u_HC45(f, i, :, :) = squeeze(u);
        [freq,PSD_HC45(f,i,:)]= getPSD(diff(pert_x_HC45(f, i, 2, 40:end))/dt,dt);
    end
end

%"Essential Tremor" - ET
pert_x_ET45 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET45 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET45 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET45 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.7; % Delay is underestimated : Delay used is 70% of the actual delay
        [x, u, x_est] = simulation(timeStab, delayError,delayError, f, x0,0.045,I);
        pert_x_ET45(f, i, :, :) = squeeze(x);
        pert_xest_ET45(f, i, :, :) = squeeze(x_est);
        pert_u_ET45(f, i, :, :) = squeeze(u);
        [freq,PSD_ET45(f,i,:)]= getPSD(diff(pert_x_ET45(f, i, 2, 40:end))/dt,dt);
    end
end


%% PSD
constantsPlots;
F = figForInkscape(19/332*86.11,10/216*64.43);

ax = subplot(2,10,1:20,'Units','centimeters');
ax.Position = [14.8,26,25.2,34.4]/10; % define your position
hold on;


% [M,I]= max(squeeze(mean(mean(PSD_ET65,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
% 
% [M,I]= max(squeeze(mean(mean(PSD_ET45,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))

[M,I]= max(squeeze(mean(mean(PSD_ET,2),1)));
% disp('The peak of the frequency is at (Hz):')
% disp(freq(I))

% Dummy plots to have the legend in black
plot45 = plot(freq,squeeze(mean(mean(PSD_HC45,2),1))/M,'Color','k','LineWidth',thickLine,'Linestyle',"--");
plot55 = plot(freq,squeeze(mean(mean(PSD_HC,2),1))/M,'Color','k','LineWidth',thickLine);
plot65 = plot(freq,squeeze(mean(mean(PSD_HC65,2),1))/M,'Color','k','LineWidth',thickLine,'Linestyle',"-.");

%Actual plots
plot(freq,squeeze(mean(mean(PSD_HC,2),1))/M,'Color',color_c,'LineWidth',thickLine);
plot(freq,squeeze(mean(mean(PSD_ET,2),1))/M,'Color',color_p,'LineWidth',thickLine)

plot(freq,squeeze(mean(mean(PSD_HC65,2),1))/M,'Color',color_c,'LineWidth',thickLine,'Linestyle',"-.");
plot(freq,squeeze(mean(mean(PSD_ET65,2),1))/M,'Color',color_p,'LineWidth',thickLine,'Linestyle',"-.")


plot(freq,squeeze(mean(mean(PSD_HC45,2),1))/M,'Color',color_c,'LineWidth',thickLine,'Linestyle',"--");
plot(freq,squeeze(mean(mean(PSD_ET45,2),1))/M,'Color',color_p,'LineWidth',thickLine,'Linestyle',"--")


xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Variation of delay duration');
hold on;
xlim([0 15]);
leg=legend([plot45, plot55, plot65], {'d=45ms', 'd=55ms', 'd=65ms'});
leg.ItemTokenSize = [30/2,18/2];
figForInkscapeSave(F,append(figurePath,'allSim_DelayDuration'))
