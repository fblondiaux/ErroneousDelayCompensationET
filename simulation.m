function [x, u, x_est, xy, L] = simulation(time_stab, delay_error1,delay_error2, pert, x0, delta, I,AB)
    % Author: Flo Blondiaux
    % Date: Jan 2024
    %    
    % Simulates one movement using the LQG framework
    %
    % time_stab = Duration of the movement [s]
    % delay_error1 = percentage of error for the delay : (1 no error, <1 :
    %   under estimated delay,). Error in the extrapolation of the state.
    % delay_error2 = error in the extrapolation of the motor commands.
    % pert = magnitude of the perturbation triggered after 10 time steps [Nm]. 
    % x0= initial state vector.
    % delta = delay for the feedback [s]
    % I = inertia [Kgm²]
    % AB = scaling factors to create mismatch in the estimation of A and B matrices
    %   and kalman gains, [1 1 1] = no errors on the estimation. 

    if nargin == 7
        %if missing define scaling factors for AB
        AB = [1 1 1];
    end
    dt = 0.005; %s
    delta_d = round(delta / dt); %discretized delay
    nbSteps = round(time_stab / dt);
    alphaNoise_c = 0.006;
    alphaNoise_d = sqrt(alphaNoise_c^2 / dt);
    noise = 1e-5;
    
    % Definition of the state
    G = 0.14; %Nms
    tau = 0.06; %ms
    Ks=0;
    Ac = zeros(8, 8);
    Ac(1:4, 1:4) = [0, 1, 0, 0;
                   -Ks/I, -G/I, 1/I, 1/I;
                   0, 0, -1/tau, 0;
                   0, 0, 0, 0];

    Bc = [0; 0; 1/tau; 0; 0; 0; 0; 0];

    % Discretization
    Ad = expm(dt * Ac);
    Bd = discreteB(Bc, Ac, dt);
    n = size(Ad, 1);

    % Cost function
    wn = [100; 0; 0; 0];
    w = wn;
    Q_n = [w(1), 0, 0, 0, -w(1), 0, 0, 0;
           0, w(2), 0, 0, 0, -w(2), 0, 0;
           0, 0, w(3), 0, 0, 0, -w(3), 0;
           0, 0, 0, 0, 0, 0, 0, 0;
          -w(1), 0, 0, 0, w(1), 0, 0, 0;
           0, -w(2), 0, 0, 0, w(2), 0, 0;
           0, 0, -w(3), 0, 0, 0, w(3), 0;
           0, 0, 0, 0, 0, 0, 0, 0];
    Q = zeros(n, n, nbSteps);

    % Stabilization time -> Qi = Qn 
    for i = 1: nbSteps
        Q(:, :, i) = Q_n;
    end
    
    %Extra constraint on the velocity for the last timestep
    w = [100; 10; 0; 0];
    Q(:,:,end)=  [w(1), 0, 0, 0, -w(1), 0, 0, 0;
           0, w(2), 0, 0, 0, -w(2), 0, 0;
           0, 0, w(3), 0, 0, 0, -w(3), 0;
           0, 0, 0, 0, 0, 0, 0, 0;
          -w(1), 0, 0, 0, w(1), 0, 0, 0;
           0, -w(2), 0, 0, 0, w(2), 0, 0;
           0, 0, -w(3), 0, 0, 0, w(3), 0;
           0, 0, 0, 0, 0, 0, 0, 0];
 
    H = eye(n); % Identity observation matrix
    R = 1e-2;

    % Backward recursion
    S = zeros(n, n, nbSteps+1);
    L = zeros(1, n, nbSteps);
    S(:, :, nbSteps+1) = Q(:, :, nbSteps);
    for k = nbSteps : -1 : 1
        L(:, :, k) = (R + Bd' * S(:, :, k + 1) * Bd)\ Bd' * S(:, :, k + 1) * Ad;
        S(:, :, k) = Q(:, :, k) + Ad' * S(:, :, k + 1) * (Ad - Bd * L(:, :, k));
    end

    % Initialization vectors
    x = zeros(n, 1, nbSteps);
    xy = zeros(n, 1, nbSteps);
    x_est = zeros(n, 1, nbSteps);
    y = zeros(n, 1, nbSteps);
    u = zeros(1, nbSteps);
    x_p = zeros(n, 1, nbSteps);
    x(:, :, 1) = x0;
    x_est(:, :, 1) = x0;

    % Noise
    C = alphaNoise_d * Bd; % Signal dependent noise
    oXi = noise * Bd * Bd';
    ceta = noise * eye(n);
    oOmega = noise * eye(n); %Sensitive noise
    Sigmaxx = noise * eye(n);

    % Forward recursion
    for k = 1 : nbSteps - 1
        if (k == 10) % Turn on perturbation (torque)
            x(4,:,k) = x(4,:,k)-pert;
        end
        % if k == 20 %Turn off perturbation
        %     x(4,:,k) = 0;
        % end

        u(:, k) = -L(:, :, k) * x_est(:, :, k); %Motor command
         
        sdn = normrnd(0, 1) * C * u(:, k); %signal dependent noise
         
        x(:, :, k + 1) = Ad * x(:, :, k) + Bd * u(:, k) + sdn + mvnrnd(zeros(n, 1), oXi)'; %True state
        y(:, :, k) = H * x(:, :, max(1,k - delta_d)) + mvnrnd(zeros(n, 1), oOmega)'; %Delayed sensory feedback
        [xy(:, :, k),V, M] = integrate(y(:, :, k), dt, delta, ...
                             u(:, max(k - delta_d, 1) : k), Ac*AB(1), Bc*AB(2), ...
                             delay_error1,delay_error2); %Extrapolation over delay interval
        x_p(:, :, k) = Ad*AB(1) * x_est(:, :, k) + Bd*AB(2) * u(:, k)+mvnrnd(zeros(n,1),ceta)'; %State prediction
        Ve = (alphaNoise_c^2 * dt * Bc * u(:, k) * u(:, k)' * Bc') + oXi + ceta;
        Sigmaxy = Ad*AB(1) * Sigmaxx;
        Sigmayy = Sigmaxx + M * oOmega * M' + alphaNoise_c^2 * V + (delta / dt) * oXi + ceta;
        Sigmaxx = Ad*AB(1) * Sigmaxx * Ad'*AB(1) + Ve + ceta - Sigmaxy/Sigmayy* Sigmaxy'; % Eq 2.13 Crevecoeur 2019
        x_est(:, :, k + 1) = x_p(:, :, k) + Sigmaxy/Sigmayy*AB(3)* (xy(:, :, k) - x_est(:, :, k)); % Eq 2.12 Crevecoeur 2019
        
    end
end
