% Script to investigate the effect of delay error on the oscillations. 
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
nbForce = 3;% length(forces);
x0 = [0;0;0;0;0;0;0;0]; %Initial state vector, state variables include 
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque

%"Healthy Controls - HC"
pert_x_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt)); % 8 = nb of state variables
pert_xest_HC = zeros(nbForce, nbSim, 8, round((timeStab) / dt));
pert_u_HC = zeros(nbForce, nbSim, 1, round(timeStab / dt)); %1 = nb of control variables
PSD_HC = zeros(nbForce,nbSim,65);

for f = 1:3
    for i = 1:nbSim
        delayError = 1; % Delay error in percentage
        [x, u, x_est,xy,L_HC] = simulation(timeStab, delayError, delayError, f, x0,delta,I);
        pert_x_HC(f, i, :, :) = squeeze(x);
        pert_xest_HC(f, i, :, :) = squeeze(x_est);
        pert_u_HC(f, i, :, :) = squeeze(u);
        [freq,PSD_HC(f,i,:)]= getPSD(diff(pert_x_HC(f, i, 2, 40:end))/dt,dt);
    end
end

% Delay error = 60%
pert_x_ET6 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET6 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET6 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET6 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.6; 
        [x, u, x_est,xy, L_ET10] = simulation(timeStab, delayError,delayError, f, x0,delta,I);
        pert_x_ET6(f, i, :, :) = squeeze(x);
        pert_xest_ET6(f, i, :, :) = squeeze(x_est);
        pert_u_ET6(f, i, :, :) = squeeze(u);
        [freq,PSD_ET6(f,i,:)]= getPSD(diff(pert_x_ET6(f, i, 2, 40:end))/dt,dt);
    end
end

%Delay error = 70%
pert_x_ET7 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET7 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET7 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET7 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.7; 
        [x, u, x_est,xy, L_ET] = simulation(timeStab, delayError,delayError, f, x0,delta,I);
        pert_x_ET7(f, i, :, :) = squeeze(x);
        pert_xest_ET7(f, i, :, :) = squeeze(x_est);
        pert_u_ET7(f, i, :, :) = squeeze(u);
        [freq,PSD_ET7(f,i,:)]= getPSD(diff(pert_x_ET7(f, i, 2, 40:end))/dt,dt);
    end
end

%Delay error = 80%
pert_x_ET8 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET8 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET8 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET8 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.8; 
        [x, u, x_est,xy, L_ET10] = simulation(timeStab, delayError,delayError, f, x0,delta,I);
        pert_x_ET8(f, i, :, :) = squeeze(x);
        pert_xest_ET8(f, i, :, :) = squeeze(x_est);
        pert_u_ET8(f, i, :, :) = squeeze(u);
        [freq,PSD_ET8(f,i,:)]= getPSD(diff(pert_x_ET8(f, i, 2, 40:end))/dt,dt);
    end
end

%Delay error = 90%
pert_x_ET9 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_ET9 = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_ET9 = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_ET9 = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=.9; 
        [x, u, x_est,xy, L_ET20] = simulation(timeStab, delayError,delayError, f, x0,delta,I);
        pert_x_ET9(f, i, :, :) = squeeze(x);
        pert_xest_ET9(f, i, :, :) = squeeze(x_est);
        pert_u_ET9(f, i, :, :) = squeeze(u);
        [freq,PSD_ET9(f,i,:)]= getPSD(diff(pert_x_ET9(f, i, 2, 40:end))/dt,dt);
    end
end




%% Power spectral density - Normalized
constantsPlots;
F = figForInkscape(19/332*86.11,10/216*64.43);
ax = subplot(2,10,1:20,'Units','centimeters');
ax.Position = [14.8,26,25.2,34.4]/10; % define your position
hold on;
% [~,I]= max(squeeze(mean(mean(PSD_ET8,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
% [~,I]= max(squeeze(mean(mean(PSD_ET9,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))
% [~,I]= max(squeeze(mean(mean(PSD_ET6,2),1)));
% disp('The peak of the frequency is at (Hz):');disp(freq(I))

[M,I]= max(squeeze(mean(mean(PSD_ET7,2),1)));
%disp('The peak of the frequency is at (Hz):')
%disp(freq(I))

plot(freq,squeeze(mean(mean(PSD_ET6,2),1))/M,'Color',color_var2,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_ET7,2),1))/M,'Color',color_p,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_ET8,2),1))/M,'Color',color_var3,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_ET9,2),1))/M,'Color',color_var1,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_HC,2),1))/M,'Color',color_c,'LineWidth',thickLine)
leg = legend(["^d = 0.6d","^d = 0.7d","^d = 0.8d","^d = 0.9d","^d = 1d"],'FontSize',4, 'Location','best');
leg.ItemTokenSize = [30/3,18/3];
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Variation of delay errors');
hold on;
xlim([0 15]);
figForInkscapeSave(F,append(figurePath,'allSim_delayErrors'))
