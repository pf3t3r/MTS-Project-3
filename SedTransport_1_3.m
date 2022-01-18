clc; clear; close all;

%% 1.3: Sensitivity of the Duration Asymmetry to ...

%% 1.3.a: Fall Velocity

alpha = 1e-4;                   % Erosion coefficent
Kv = 1e-2;                      % Vertical eddy diffusivity (vertical mixing)
WS = linspace(0.5e-3,2e-2,5);   % Possible values for the fall velocity
diffWs = [];

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
    
    % Save C and U for each value of WS at the first point in x.    
    C_Ws(i,:) = C(1,:);
    U_Ws(i,:) = U(1,:);
    
    firstPeakInC = max(C(1,100:200));
    firstPeakInU = max(U(1,100:200));
    
    timeOfFirstPeakInC = t(C(1,:)==firstPeakInC);
    timeOfFirstPeakInU = t(U(1,:)==firstPeakInU);
    
    durationDiff = (timeOfFirstPeakInC - timeOfFirstPeakInU)/3600;
    diffWs = [diffWs durationDiff];
end

% Quick update for legends
 for i = 1:length(WS)
      Ws_legend{i} = num2str(WS(i),'W_s = %.4f m/s');
 end

figure

subplot(4,1,1)
yyaxis left
plot(t,C_Ws(1,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Ws(1,:))
ylabel('U [m/s]');
% legend('Ws = 0.0005');
legend(Ws_legend(1));
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Fall Velocity: I')

subplot(4,1,2)
yyaxis left
plot(t,C_Ws(2:3,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Ws(2:3,:))
ylabel('U [m/s]');
legend(Ws_legend(2:3));
grid(gca,'minor');
grid on;
title('Sensitivity of Duration Asymmetry to Fall Velocity: II - III')

subplot(4,1,3)
yyaxis left
plot(t,C_Ws(4:end,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Ws(4:end,:))
ylabel('U [m/s]');
legend(Ws_legend(4:end));
grid(gca,'minor');
grid on;
title('Sensitivity of Duration Asymmetry to Fall Velocity: IV - X')

subplot(4,1,4)
plot(WS,diffWs)
xlabel('W_s [m/s]');
ylabel('\Deltat [hrs]');
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Fall Velocity')

savefig('Matlab3_1_3_i');


%% 1.3.b: Eddy Diffusivity

alpha = 1e-4;                   % Erosion coefficent
Ws = 1e-3;                      % Fall velocity of sediment
KV = linspace(10e-3,10e-1,5);    % Array of eddy diffusivities
diffKv = [];

% Iterate across possible eddy diffusivities
for i = 1:length(KV)
    
    Kv = KV(i);
    
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
    
    % Save C and U for each value of KV at the first point in x.    
    C_Kv(i,:) = C(1,:);
    U_Kv(i,:) = U(1,:);
    
    firstPeakInC = max(C(1,100:200));
    firstPeakInU = max(U(1,100:200));
    
    timeOfFirstPeakInC = t(C(1,:)==firstPeakInC);
    timeOfFirstPeakInU = t(U(1,:)==firstPeakInU);
    
    durationDiff = (timeOfFirstPeakInC - timeOfFirstPeakInU)/3600;
    diffKv = [diffKv durationDiff];
  
end

% Quick update for legends
 for i = 1:length(KV)
      Kv_legend{i} = num2str(KV(i),'K_v = %.4f m^2/s');
 end

figure

subplot(4,1,1)
yyaxis left
plot(t,C_Kv(1,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Kv(1,:))
ylabel('U [m/s]');
legend(Kv_legend(1));
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Eddy Diffusivity: I - III')

subplot(4,1,2)
yyaxis left
plot(t,C_Kv(2:3,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Kv(2:3,:))
ylabel('U [m/s]');
legend(Kv_legend(2:3));
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Eddy Diffusivity: IV - VI')

subplot(4,1,3)
yyaxis left
plot(t,C_Kv(4:end,:))
xlabel('Time [s]');
ylabel('C [kg/m^2]');
hold on
yyaxis right
plot(t,U_Kv(4:end,:))
ylabel('U [m/s]');
legend(Kv_legend(4:end));
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Eddy Diffusivity: VII - X')

subplot(4,1,4)
plot(KV,diffKv)
xlabel('K_v [m^{2}/s]');
ylabel('\Deltat [hrs]');
grid(gca,'minor')
grid on;
title('Sensitivity of Duration Asymmetry to Eddy Diffusivity')

savefig('Matlab3_1_3_ii');