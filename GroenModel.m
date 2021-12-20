function [Conc]=Groenmodel(U,t,deltaT,T,Ws,alpha,Kv)

Tend=t(end);
deltaTfix=Kv/Ws^2;            % Maximum time step allowed.

if deltaTfix>deltaT             % If deltaT<deltaTmax 
    deltaTfix=deltaT;           % WHY?
else
    deltaT=deltaTfix;
end

% Within the code deltaT can be time dependent. Most often deltaTfix will
% be used. 

% initialise solution at t=0
k=1;
tt(1)=t(1);
C(1)=0.0; %concentration is zero at the start of calculations 

while tt(k)<Tend  
    %************************************************
    
    % Predictor step
    Uf(k) = interp1(t,U,tt(k));           % Because time step might be different than prescribed, an interpolation step is needed.
   
    % Make rest of code here. See assignment for more information.
    deltaT(k) = deltaT;
    E(k) = alpha * U.^2;
    D(k) = Ws.^2 ./ Kv;
    % New C will be caused by difference betwen erosion and deposition
    C(k+1) = C(k) + (E(k) - D(k))*deltaT(k);
    
    %*********************************************
    % End of predictor step
    %*********************************************
    
    % Problem: predicted sediment concentration can become negative. Check for that. If C is negative: take smaller
    % time step.
    if C(k+1)<0
        deltaT(k)=C(k)/D(k)*0.5;            % C(k)/D(k) is estimated time to deposit all sediment that is in water column. 
        C(k+1)=C(k)+(E(k)-D(k))*deltaT(k);  % 
    else
        deltaT(k)=deltaTfix;
    end
    
    if tt(k)+deltaT(k)>Tend                 % Don't let the time go beyond the final time step. 
        deltaT(k)=Tend-tt(k);
    end
   
    % Go one step further. Update time and k.
    tt(k+1)=tt(k)+deltaT(k);
    
    %*****************************************************************
    % corrector step. see assignment. 
    
    % Uf(k+1)=interp1
    % E(k+1)=
    % Etc
    
    
    
    %******************************************************************
    % End of corrector step
    
    k=k+1;
    
    
end

Conc=interp1(tt,C,t);           % Since we have calculated the solution on a new time vector, we have to interpolate the solution to the right time vector. 


