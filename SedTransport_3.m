clc; clear; close all;

%% Analysis of Sediment Transport in the Gironde Estuary

% Coefficients and typical values
alpha = 1e-4;                   % Erosion coefficent
Ws = 3e-3;                      % Fall velocity
Cd = 0.2e-3;                    % Typical drag coefficient [] (Matlab 2)
U = 8;                          % Typical flow velocity [m/s] (Matlab 2)
depth = 8.5;                    % Equilibrium Depth [m] (Matlab 2) 
Kv = Cd*U*depth;                % Vertical Eddy Diffusivity; this should be 
                                % close to 1e-2.

% Gironde Estuary: spatial parameters [width, length, etc.]
xr = 110e3;              % Point at which estuary stops narrowing [m]
                         % In our case, the estuary splits into two rivers.
                         % The width becomes constant after 107 and 111 km
                         % respectively, so we pick 110 km as an average.
Lbasin = 1.5*xr;         % Length of the basin / estuary [m]
B0 = 5.8e3;              % Width of the basin at seaward side [m]
                         % Based on Google Earth view of the Gironde mouth.
Briver = 0.25e3;         % Width of the river [m] where estuary stops
                         % narrowing: 0.2 km on Dordogne, 0.3 km on Garonne,
                         % => average taken. 
deltaX = Lbasin/80;      % Spatial step [m]
Lb = -xr/log(Briver/B0); % E-folding length scale [m]   
M2amp = 1.55;            % Amplitude of M2 tide at seaward side [m]
                         % This is the water level at the mouth.
discharge = 0;           % Constant river discharge at landward boundary. 
                         % Prescribed as zero for simplicity in the project
                         % description.
x=0:deltaX:Lbasin;       % With deltax = c.1.3e3, we have 80 grid points.

% Water levels forced with D2 tide only
Td2 = 12*3600 + 25*60;   % M2 tidal period [s]
deltaT = 25;             % Time step [s]. Must satisfy Courant condition.
time = 0:deltaT:30*Td2;  % Time [s]
Nt = length(time);       % Size of time array

% Define frequencies to be analysed. M1 tide is not included.
global wn
wn(1)=2*pi/Td2;         % M2
wn(2)=2*wn(1);          % M4
wn(3)=3*wn(1);          % M6


% Evaluate velocity asymmetry for different phases
% PhaseM4=(0:45:180)/180*pi;
% for i = 1:length(PhaseM4)
    
% phaseM4=PhaseM4(i);
% Maybe this is not correct.
phaseM4 = 60/180*pi;      % Phase = 60 deg.

T = (12*60+25)*60;        % M2 and M4 tide. Time in seconds. 
Tend = 10*T;              % Five tidal periods modeled -> for very fine
                          % sand and large erosion constants more tidal
                          % periods need to be solved
deltaT = 300;             % Time step of 5 minutes
t = 0:deltaT:Tend;
Nt = length(t);

ampD1 = 0;              % D1 does not need to be prescribed here.
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

% L = 1e4;                 % We model a simple basin with a length of ten km
L = Lbasin;
% dx = 400;                % Grid distance
dx = deltaX;
% x = 0:dx:L;              % x-coordinate. Seaward end is at x = L,
                         % landward end at x=0. 
Nx = length(x);                       

% H = 10-8e-4*x;           % Bottom profile. Linear sloping bottom. 2 m deep
                         % near landward boundary, 10 m deep near inlet. 
H = depth*ones(1,length(x));
dHdx(1:Nx) = -8e-4;

%**************************************************************************
% After a call to hydromodel, flow velocity at each position as a
% function of time is known
%**************************************************************************

U = HydroModel2(t, Z, dZdt, H, dHdx, x, dx);

%**************************************************************************
% Here the sediment concentrations are calculated for each position in the
% basin with the Groen model.

for px = 1:Nx
    [C(px, 1:Nt)] = GroenModel(U(px, 1:Nt), t, deltaT, T, Ws, alpha, Kv);
end

% Save the value of C and U at the first point in X    
%     C_phaseM4(i,:) = C(1,:);
C_x1 = C(1,:);
%     U_phaseM4(i,:) = U(1,:);
U_x1 = U(1,:);

% Determine tidally-averaged sediment transport as a function of position
% also U coz y not m8
for px = 1:Nx
    C_avg(px) = mean(C(px,:));
    U_avg(px) = mean(U(px,:));
end
disp(C_avg)

figure

subplot(3,1,1)
yyaxis left
plot(t/3600,C_x1);
ylabel('C [kg/m^2]');
yyaxis right
plot(t/3600,U_x1);
hold off
ylabel('U [m/s]');
xlabel('t [hrs]');
% legend(PhaseM4_legend);
grid(gca,'minor')
grid on;
title('Velocity and Sediment Concentration at the first location in X')

subplot(3,1,2)
yyaxis left
plot(t(1:300)/3600,C_x1(:,1:300));
ylabel('C [kg/m^2]');
yyaxis right
plot(t(1:300)/3600,U_x1(:,1:300));
ylabel('U [m/s]');
xlabel('t [hrs]');
% legend(PhaseM4_legend);
grid(gca,'minor')
grid on;
title('First 300 values of U and C at the first point in X')

subplot(3,1,3)
yyaxis left
plot(x,U_avg);
ylabel('U [m/s]');
hold on
yyaxis right
plot(x,C_avg);
ylabel('C [kg/m^2]');
hold off
title('C and U as a function of X');

savefig('pt-3-i');