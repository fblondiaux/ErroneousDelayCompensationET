function Bd = discreteB(B, A, dt)
    % The routine uses a dt/100 discretization step
    % Correspond to eq. 9 of the publication
    step = dt/100;
    delta = 0:step:dt; 
    Bd = 0;
    for i = 1:size(delta,2)
        Bd = Bd + step* expm(delta(i) * A);
    end
    Bd = Bd * B;
end
