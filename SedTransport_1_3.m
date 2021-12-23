clc; clear; close all;

%% 1.3: Sensitivity of the Duration Asymmetry to ...

%% 1.3.a: Fall Velocity

alpha = 1e-4;                   % Erosion coefficent
Kv = 1e-2;                      % Vertical eddy diffusivity (vertical mixing)
WS = linspace(0.5e-3,2e-2,5);   % Possible values for the fall velocity
 
durationAsymmetry = [];         % Difference between peak ebb and flood

% Iterate across possible fall velocities

for i = 1:length(WS)
    
    Ws = WS(i);
    
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
    phaseM4 = pi/2;

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
    
    test(i,:) = C(1,:);
    test2(i,:) = U(1,:);
    
    % Peak sediment concentration at flood and ebb
    C_Max = findpeaks(C(1,201:371));
    
    % Minimum sediment concentration at flood and ebb
    C_MinTF = islocalmin(C(1,201:371));
    C_Min = nonzeros(C_MinTF.*C(1,201:371));
    % Time at peak flood and ebb
    
    % Difference between peaks
    if length(C_Max)==1
        Diff = abs(C_Max(1) - max(C(1,268:371)));
    else
        Diff = C_Max(1) - C_Max(2);
    end
    durationAsymmetry = [durationAsymmetry Diff];
end


%Tracking the behaviours and locating the peaks
figure
yyaxis left
plot(t,test(2:5,:))
ylabel('C [..]');
% legend('Ws = 0.0103','Ws = 0.0151','Ws = 0.0200');
hold on
yyaxis right
plot(t,test2(2:5,:))
ylabel('U [m/s]');
grid on;


figure
plot(t,U);
grid on;

% figure
% scatter(WS,durationAsymmetry);
% grid on;

% hold on
% for i = 1:length(peaks)-1
%     xline(peaks(i),'lineWidth',10);
% end
% hold off
% plot(WS,abs(DIFF_ws_3))
% title('Sensitivity analysis of duration asymmetry to fall velocity (W_s)');
% xlabel('W_{s} [m/s]');
% ylabel('Concentration [kg/m^2]');
% grid on;
% savefig('Matlab3_1_iv');


%% 1.3.b: Eddy Diffusivity

alpha = 1e-4;                   % Erosion coefficent
Ws = 1e-3;                % Fall velocity of sediment
KV = linspace(1e-3,1e-1,5);