function [Conc]=Cmodel(C,U,t,deltaT, T, Ws, alpha, Kv, dQsdx)

Tend=t(end);
deltaTfix=Kv/Ws^2;            % Maximum time step allowed.

if deltaTfix>deltaT             % If deltaT<deltaTmax 
    deltaTfix=deltaT;           
else
    deltaT=deltaTfix;
end

% Within the code deltaT can be time dependent. Most often deltaTfix will
% be used. 

% initialise solution at t=0
k=1;
tt(1)=t(1);
C2(1)=C(1); %concentration is equal to concentration from GroenModel at the start of calulations
A(1) = dQsdx(1); %advection at start of calculations from GroenModel

while tt(k)<Tend  
    %************************************************
    
    % Predictor step
    Uf(k) = interp1(t,U,tt(k));           % Because time step might be different than prescribed, an interpolation step is needed.
   
    % Define variables
    deltaT(k) = deltaTfix;
    E(k) = alpha * Uf(k).^2;
    D(k) = (Ws.^2 ./ Kv) * C2(k);
    Qs(k) = Uf(k).*C2(k);
    dx = 400;
    A(k) = gradient(Qs(k),deltaT,dx);
    
    % New C will be caused by difference betwen erosion and deposition with
    % advection. It includes a d/dx UC term
    
    C2(k+1) = C2(k) + (E(k) - D(k) - A(k))*deltaT(k);
    
    %*********************************************
    % End of predictor step
    %*********************************************
    
    % Problem: predicted sediment concentration can become negative. Check
    % for that. If C is negative: take smaller time step.
    if C2(k+1)<0
        deltaT(k)=C2(k)/D(k)*0.5;            % C(k)/D(k) is estimated time
                                            % to deposit all sediment that
                                            % is in water column. 
        C2(k+1)=C2(k)+(-A(k)+E(k)-D(k))*deltaT(k);  % 
    else
        deltaT(k)=deltaTfix;
    end
    
    if tt(k)+deltaT(k)>Tend                 % Don't let the time go beyond
                                            % the final time step. 
        deltaT(k)=Tend-tt(k);
    end
   
    % Go one step further. Update time and k.
    tt(k+1)=tt(k)+deltaT(k);
    
    %*****************************************************************
    % Corrector step. 
    
    % Define variables 
    Uf(k+1)=interp1(t,U,tt(k+1));
    E(k+1) = alpha * Uf(k+1).^2;
    D(k+1) = (Ws.^2./Kv) * C2(k+1);
    Qs(k+1) = Uf(k+1).*C2(k+1);
    A(k+1) = gradient(Qs(k+1),deltaT,dx);
    
    % New C will be caused by difference betwen erosion and deposition
    C2(k+1) = C2(k) + (E(k) + E(k+1) - D(k) - D(k+1) - A(k) - A(k+1)) * (deltaT(k)*0.5);
    
    %******************************************************************
    % End of corrector step
    
    k=k+1;
    
end

Conc=interp1(tt,C2,t);           % Since we have calculated the solution on a new time vector, we have to interpolate the solution to the right time vector. 


