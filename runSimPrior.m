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
nbSim = 20;
timeStab = 0.8; %Stabilization time
dt = 0.005; %[s]
delta = 0.055; %Delay of the feedback [s]
I = 0.15; 
forces = [2]; % [Nm]
nbForce = 3;% length(forces);
x0 = [0;0;0;0;0;0;0;0]; %Initial state vector, state variables include 
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque

%"Healthy Controls - HC"
pert_x_HC = zeros(nbForce, nbSim, 8, round(timeStab / dt)); % 8 = nb of state variables
pert_xest_HC = zeros(nbForce, nbSim, 8, round(timeStab / dt));
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

% under - Error with un underestimation of the prior = 0.7prior
pert_x_underEst = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_underEst = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_underEst = zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_underEst = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=1; % Delay is underestimated : Delay used is 70% of the actual delay
        [x, u, x_est,xy, L_ET] = simulation(timeStab, delayError,delayError, f, x0,delta,I,[0.7 1]);
        pert_x_underEst(f, i, :, :) = squeeze(x);
        pert_xest_underEst(f, i, :, :) = squeeze(x_est);
        pert_u_underEst(f, i, :, :) = squeeze(u);
        [freq,PSD_underEst(f,i,:)]= getPSD(diff(pert_x_underEst(f, i, 2, 40:end))/dt,dt);
    end
end

% B - Error only in B with Best = 1/2B
pert_x_overEst = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_xest_overEst = zeros(nbForce, nbSim, 8, round(timeStab / dt));
pert_u_overEst= zeros(nbForce, nbSim, 1, round(timeStab / dt));
PSD_overEst = zeros(nbForce,nbSim,65);
%With a error
for f = 1:3
    for i = 1:nbSim
        delayError=1; % Delay is underestimated : Delay used is 70% of the actual delay
        [x, u, x_est,xy, L_ET] = simulation(timeStab,delayError,delayError,f, x0,delta,I,[1.1 1]);
        pert_x_overEst(f, i, :, :) = squeeze(x);
        pert_xest_overEst(f, i, :, :) = squeeze(x_est);
        pert_u_overEst(f, i, :, :) = squeeze(u);
        [freq,PSD_overEst(f,i,:)]= getPSD(diff(pert_x_overEst(f, i, 2, 40:end))/dt,dt);
    end
end


%% Plot with the arm angle, angular velocity, control and PSD

tv = -10*dt:dt:(timeStab-dt-10*dt); % Perturbation start after 10 timesteps
tv=tv*1000; %From s to ms

constantsPlots;

%Define the size of the figure
F = figForInkscape(19/332*86.11,11/216*64.43);
ax = subplot(1,11,1:5,'Units','centimeters');
ax.Position = [14.8,66.5,33.5,23.86]/10; 
hold on;
for f = forces
    plot(tv,squeeze(mean(pert_x_HC(f, :, 1, :), 2))*180/pi,'Color',color_c,'LineWidth',thickLine);
    plot(tv,squeeze(mean(pert_x_underEst(f, :, 1, :), 2))*180/pi,'Color',color_var1,'LineWidth',thickLine); %,
    plot(tv,squeeze(mean(pert_x_overEst(f, :, 1, :), 2))*180/pi,'Color',color_var2,'LineWidth',thickLine); %   
end
xline(0);
xlabel('Time (ms)');
ylabel('Angle (deg)');
title('Elbow joint angle');
xlim([tv(1) tv(end)])

%For different error of the Kalman gains
%leg = legend(["No Error","K=0.2*K","K=2.5*K",],'FontSize',4, 'Location','best');

leg = legend(["No Error","prior=0.7*prior","prior=1.1*prior",],'FontSize',4, 'Location','best');
leg.ItemTokenSize = [30/3,18/3];
%figForInkscapeSave(F,append(figurePath,'runSimAngle'))


%% Power spectral density - Normalized
ax = subplot(1,11,6:7,'Units','centimeters');
ax.Position = [56.3,66.5,15.2,23.86]/10; % define your position
hold on;

[M,I]= max(squeeze(mean(mean(PSD_overEst(:,:,:),2),1)));
disp('The peak of the frequency is at (Hz):')
disp(freq(I))

[M,I]= max(squeeze(mean(mean(PSD_underEst(:,:,:),2),1)));
disp('The peak of the frequency is at (Hz):')
disp(freq(I))



plot(freq,squeeze(mean(mean(PSD_HC,2),1))/M,'Color',color_c,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_underEst,2),1))/M,'Color',color_var1,'LineWidth',thickLine)
plot(freq,squeeze(mean(mean(PSD_overEst,2),1))/M,'Color',color_var2,'LineWidth',thickLine)
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
ax = subplot(1,11,8:9,'Units','centimeters');
ax.Position = [79.5,66.5,13,23.86]/10; % define your position
hold on;
errorbar(forces,squeeze(mean(mean(pert_u_HC(:,:,1,19:25),4),2))',[0 0 0],'-o','Color',color_c,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor', color_c,'MarkerSize',5);
errorbar(forces,squeeze(mean(mean(pert_u_underEst(:,:,1,19:25),4),2))',[0 0 0],'-o','Color',color_var1,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor',color_var1,'MarkerSize',5);
errorbar(forces,squeeze(mean(mean(pert_u_overEst(:,:,1,19:25),4),2))',[0 0 0],'-o','Color', color_var2,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor',color_var2,'MarkerSize',5);
axis([0.3 3.5 -5 25])
xlabel('Perturbation (Nm)')
ylabel('Control (a.u.)')
title('R2')

%LLR - R3
ax = subplot(1,11,10:11,'Units','centimeters');
ax.Position = [100.5,66.5,13,23.86]/10; % define your position
hold on;
errorbar(forces,squeeze(mean(mean(pert_u_HC(:,:,1,25:31),4),2))',[0 0 0],'-o','Color',color_c,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor', color_c,'MarkerSize',5);
errorbar(forces,squeeze(mean(mean(pert_u_underEst(:,:,1,25:31),4),2))',[0 0 0],'-o','Color',color_var1,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor', color_var1,'MarkerSize',5);
errorbar(forces,squeeze(mean(mean(pert_u_overEst(:,:,1,25:31),4),2))',[0 0 0],'-o','Color', color_var2,'LineWidth',thickLine,'CapSize',0,'MarkerFaceColor', color_var2,'MarkerSize',5);
axis([0.5 3.5 -5 25])
title('R3')

sgtitle('Impact of an error in the internal prediction')
figForInkscapeSave(F,append(figurePath,'allSim_prior'))
