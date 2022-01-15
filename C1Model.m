function [Conc] = CModel(U, t,deltaT, T, Ws, alpha, Kv, dQsdx)
Tend = t(end);
deltaTfix = Kv/Ws^2; %maximum time step allowed.

%Within the code deltaT can be time dependent. Most often deltaTfix will be used.
if deltaTfix > deltaT %if deltaT < deltaTmax 
    deltaTfix = deltaT;
else
    deltaT = deltaTfix;
end

%Initialise solution at t=0
k = 1;
tt(1) = t(1);
C(1) = 0.0;

while tt(k) < Tend  
    %Predictor step
    Uf(k) = interp1(t, U, tt(k)); %because time step might be different than prescribed, an interpolation step is needed.
    dQsdxf(k) = interp1(t, dQsdx, tt(k)); %because time step might be different than prescribed, an interpolation step is needed.
    deltaT(k) = deltaTfix;
    E(k) = alpha * Uf(k)^2;
    D(k) = (Ws^2/Kv) * C(k);
    C(k+1) = C(k) + (E(k)-D(k) - dQsdxf(k)) * deltaT(k); %new C will be caused by difference betwen erosion and deposition
        
    %Problem: predicted sediment concentration can become negative. Check for that. If C is negative: take smaller time step.
    if C(k+1) < 0
        deltaT(k) = C(k)/D(k) * 0.5; %C(k)/D(k) is estimated time to deposit all sediment that is in water column. 
        C(k+1) = C(k) + (E(k)-D(k) - dQsdxf(k)) * deltaT(k);
    else
        deltaT(k) = deltaTfix;
    end
    
    if tt(k) + deltaT(k) > Tend %don't let the time go beyond the final time step. 
        deltaT(k) = Tend - tt(k);
    end
   
    %Go one step further. Update time and k.
    tt(k+1) = tt(k) + deltaT(k);
    
    %Corrector step
    %Uf(k+1) = interp1(t, U, tt(k+1));
    %dQsdxf(k) = interp1(t, dQsdx, tt(k));
    %E(k+1) = alpha * Uf(k)^2 + alpha * Uf(k+1)^2;
    %D(k+1) = (Ws^2/Kv) * C(k) + (Ws^2/Kv) * C(k+1);
    %C(k+1) = C(k) + (E(k+1)-D(k+1)) * 0.5*deltaT(k);

    k = k+1;
end

Conc = interp1(tt,C,t); %since we have calculated the solution on a new time vector, we have to interpolate the solution to the right time vector. end