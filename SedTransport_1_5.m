clc; clear; close all;

%% 1.4: Sensitivity of tidally averaged sediment transport to the relative phase difference between M2 and M4

alpha = 1e-4;                   % Erosion coefficent
Kv = 10e-2;                     % Vertical eddy diffusivity (vertical mixing)
Ws = 1e-3;                      % Possible values for the fall velocity

% Evaluate velocity asymmetry for different phases
PhaseM4=(0:45:180)/180*pi;
for i = 1:length(PhaseM4)
    
    phaseM4=PhaseM4(i);
    
    T = (12*60+25)*60;        % M2 and M4 tide. Time in seconds. 
    Tend = 10*T;              % Five tidal periods modeled -> for very fine
                              % sand and large erosion constants more tidal
                              % periods need to be solved
    deltaT = 300;             % Time step of 5 minutes
    t = 0:deltaT:Tend;
    Nt = length(t);
    
    ampD1 = 0;                % in part 1 and 2 D1 = 0. Depending on your
                              % estuary, you might want to prescribe D1 for
                              % part 3. 
    ampM2 = 1;
    ampM4 = 0.2;
    phaseD1 = 0;
    phaseM2 = 0;

    % Water level prescribed below as a sine function.
    Z = ampD1*sin(pi*t/T + phaseD1) + ampM2*sin(2*pi*t/T + phaseM2) + ...
        ampM4*sin(4*pi*t/T + phaseM4);   
    
    % Flow velocity will behave as a cosine function.
    dZdt = ampD1*1*pi/T*cos(pi*t/T + phaseD1) + ... 
           ampM2*2*pi/T*cos(2*pi*t/T + phaseM2) + ...
           ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); 

    L = 1e4;                 % We model a simple basin with a length of ten km
    dx = 400;                % Grid distance
    x = 0:dx:L;              % x-coordinate. Seaward end is at x=L, landward end at x=0. 
    Nx = length(x);                       

    H = 10-8e-4*x;           % Bottom profile. Linear sloping bottom. 2 m deep
                             % near landward boundary, 10 m deep near inlet. 
    dHdx(1:Nx) = -8e-4;

    %**************************************************************************
    % After a call to hydromodel, flow velocity at each position as a
    % function of time is known
    %**************************************************************************

    U = HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

    %**************************************************************************
    % Here the sediment concentrations are calculated for each position in the
    % basin with the Groen model.

    for px=1:Nx
        [C(px,1:Nt)] = GroenModel(U(px,1:Nt), t, deltaT, T, Ws, alpha, Kv);
    end
    
    Qs=U.*C;                                            % Qs is sediment flux

    Nsteps=T/deltaT;                                    % Nr of timestepf in one tidal cycle.

    % Calculate tidally averaged sediment transport as a function of position
    % in the estuary (only averaging over last tidal cycle)
    Qs_x = [];
    U_x = [];
for position = 1:26
    Qs_t(position)=0;
    U_t(position)=0;
for time = 1342:1491
    Qs_t(position) = Qs_t(position) + Qs(position,time);
    U_t(position) = U_t(position) + U(position,time);
end
    % Tidally-averaged sediment transport
    meanQs_x = Qs_t(position)/149;   
    Qs_x = [Qs_x meanQs_x];

    meanU_x = U_t(position)/149;
    U_x = [U_x meanU_x];
end      
    Qs_X(i,:)=Qs_x;
    U_X(i,:)=U_x;
end

% Quick update for legends
 for i = 1:length(PhaseM4)
      PhaseM4_legend{i} = num2str(PhaseM4(i)/pi*180,'phase = %.0f°');
 end

figure
plot(x/1000,Qs_X*1000);
ylabel('Q_s [g m^{-1} s^{-1}]');
xlabel('x [km]');
legend(PhaseM4_legend);
grid(gca,'minor')
grid on;
title('Sensitivity of tidally averaged sediment transport to the relative phase difference between M2 and M4');
savefig('Matlab3_1_5_i');