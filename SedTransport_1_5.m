clc; clear; close all;

%% 1.4: Velocity Asymmetry

alpha = 1e-4;                   % Erosion coefficent
Kv = 1e-2;                      % Vertical eddy diffusivity (vertical mixing)
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
    
    % Save C and U for each value of phaseM4 at the first point in x.    
    C_phaseM4(i,:) = C(1,:);
    U_phaseM4(i,:) = U(1,:);
    
    peakFlood(i) = max(findpeaks(U_phaseM4(i,:)));
    peakEbbLoc = islocalmin(U_phaseM4(i,:));
    peakEbb(i) = -min(nonzeros(peakEbbLoc.*U_phaseM4(i,:)));
    
    netVelocity(i) = peakFlood(i) - peakEbb(i);
end

% Quick update for legends
 for i = 1:length(PhaseM4)
      PhaseM4_legend{i} = num2str(PhaseM4(i)/pi,'phase = %.2f rad');
 end

figure
subplot(3,1,1)

yyaxis left
plot(t/3600,C_phaseM4);
ylabel('C [kg/m^2]');
yyaxis right
plot(t/3600,U_phaseM4);
hold off
ylabel('U [m/s]');
xlabel('t [hrs]');
legend(PhaseM4_legend);
grid(gca,'minor')
grid on;
title('Velocity asymmetry and sediment concentration I')

subplot(3,1,2)
yyaxis left
plot(t(1:300)/3600,C_phaseM4(:,1:300));
ylabel('C [kg/m^2]');
yyaxis right
plot(t(1:300)/3600,U_phaseM4(:,1:300));
ylabel('U [m/s]');
xlabel('t [hrs]');
legend(PhaseM4_legend);
grid(gca,'minor')
grid on;
title('Velocity asymmetry and sediment concentration II')

subplot(3,1,3)
plot(PhaseM4/pi,netVelocity);
ylabel('Peak Flood - Peak Ebb [m/s]');
xlabel('Phase of M4 [\pi]');
title('Velocity asymmetry and phase');
grid on;

% A positive value for netVelocity indicates a net transport in the flood
% direction whereas a negative value indicates a net transport in the ebb
% direction.

savefig('pt-1-4');


% This figure is not needed for the report.

figure
subplot(2,1,1);
plot(t,U);
xlabel('time [s]');
ylabel('U [m/s]');
title('Velocity vs Time for 26 locations across the basin');
grid on;
legend();

subplot(2,1,2);
plot(t(104:255),U(:,104:255));
xlabel('time [s]');
ylabel('U [m/s]');
grid on;
legend();

% The highest range of velocities in the above are for values of X close
% to the seaward end of the basin. At the landward end of the basin, the
% velocities reach zero.
% The asymmetry of ebb and flood velocities stays constant across the
% basin. The overall range simply declines.