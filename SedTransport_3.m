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

T = (12*60+25)*60;        % M2 and M4 tide. Time in seconds. 
Tend = 10*T;              % Five tidal periods modeled -> for very fine
                          % sand and large erosion constants more tidal
                          % periods need to be solved
deltaT = 300;             % Time step of 5 minutes
t = 0:deltaT:Tend;
Nt = length(t);

ampD1 = 0;                % D1 does not need to be prescribed here.
ampM2 = 1;
ampM4 = 0.2;
phaseD1 = 0;
phaseM2 = 0;
phaseM4 = 60/180*pi;      % Phase = 60 deg. Might not be correct

% Water level prescribed below as a sine function.
Z = ampD1*sin(pi*t/T + phaseD1) + ampM2*sin(2*pi*t/T + phaseM2) + ...
    ampM4*sin(4*pi*t/T + phaseM4);   

% Flow velocity will behave as a cosine function.
dZdt = ampD1*1*pi/T*cos(pi*t/T + phaseD1) + ... 
       ampM2*2*pi/T*cos(2*pi*t/T + phaseM2) + ...
       ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); 

L = Lbasin;                % Length of Gironde Estuary
dx = deltaX;               % Grid spacing
Nx = length(x);                       

% Not sure if the following is correct.
% H = 10-8e-4*x;         % Bottom profile. Linear sloping bottom. 2 m deep
                         % near landward boundary, 10 m deep near inlet. 
H = depth*ones(1,length(x));
dHdx(1:Nx) = -8e-4;

%**************************************************************************
% Call hydromodel in order to find the flow velocity at each position as a
% function of time
%**************************************************************************

U = HydroModel2(t, Z, dZdt, H, dHdx, x, dx);

%**************************************************************************
% Then calculate the sediment concentrations for each position in the
% basin with the Groen model.
%**************************************************************************

for px = 1:Nx
    [C(px, 1:Nt)] = GroenModel(U(px, 1:Nt), t, deltaT, T, Ws, alpha, Kv);
end

% Save the value of C and U at the first point in X.
% This code may not be necessary.
C_x1 = C(1,:);
U_x1 = U(1,:);

% Try the same as above, but for each point in X.
% No point in doing this since we don't change the array at all.
% for px = 1:Nx
%     C_x(px) = C(px,:);
%     U_x(px) = U(px,:);
% end

% Determine tidally-averaged sediment transport as a function of position
% (and also tidally-averaged velocity)
for px = 1:Nx
    C_avg(px) = mean(C(px,:));
    U_avg(px) = mean(U(px,:));
end
disp(C_avg)

 for i = 1:Nx
      C_legend{i} = num2str(x(i)/1000,'x = %.0f km');
      U_legend{i} = num2str(x(i)/1000,'x = %.0f km');
 end

figure

subplot(3,1,1)
yyaxis left
for i = 1:int16(Nx/10):Nx
    plot(t/3600,C(i,:));
    hold on
end
ylabel('C [kg/m^2]');
legend(C_legend(1:int16(Nx/10):Nx));
yyaxis right
for i = 1:int16(Nx/10):Nx
    plot(t/3600,U(i,:));
end
legend(U_legend(1:int16(Nx/10):Nx));
hold off
ylabel('U [m/s]');
xlabel('t [hrs]');
grid(gca,'minor')
grid on;
title('Velocity and Sediment Concentration as a function of time')

subplot(3,1,2)
yyaxis left
plot(t(1:300)/3600,C_x1(:,1:300));
ylabel('C [kg/m^2]');
yyaxis right
plot(t(1:300)/3600,U_x1(:,1:300));
ylabel('U [m/s]');
xlabel('t [hrs]');
grid(gca,'minor')
grid on;
legend('C','U');
% title('U and C as a function of time: first two tidal cycles, at x = %d',x(1));
title(['U and C as a function of time at x = ',num2str(x(1)),' km, first two tidal cycles'])

subplot(3,1,3)
yyaxis left
plot(x,U_avg);
ylabel('U [m/s]');
hold on
yyaxis right
plot(x,C_avg);
ylabel('C [kg/m^2]');
hold off
grid(gca,'minor')
grid on;
legend('U_{avg}','C_{avg}');
title('Tidally-averaged Sediment Concentration and Velocity VS. basin length');

savefig('pt-3-i');