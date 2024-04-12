%% Simulation parameters
nbSim = 20;
nbState = 8; % nb of state variables
nbControl = 1; % nb of control variables
timeStab = 0.8; %Stabilization time
dt = 0.005; %[s] discretization
delta = 0.055; %Delay of the feedback [s]
delayError1 = 1; % Delay error == 100%
delayError2 = 0.7; % Delay error == 70%
I = 0.15; %Inertia [Kgm²]
forces = [1 2 3]; % [Nm]
nbForce = length(forces);
x0 = [0; 0; 0; 0; 0; 0; 0; 0]; %Initial state vector, state variables:
% angle, angular velocity, torque, external torque, target angle, target
% angular velocity, target torque and dummy variable target external torque