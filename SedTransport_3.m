clc; clear; close all;

%% Analysis of Sediment Transport in the Gironde Estuary

% Coefficients and typical values
alpha = 1e-4;                   % Erosion coefficent
Ws = 3e-3;                      % Fall velocity
Cd = 2.5e-2;                    % Typical drag coefficient [] (Matlab 2)
U = 0.052;                      % Typical flow velocity [m/s] (Matlab 2)
H0 = 8.5;                       % Equilibrium Depth [m] (Matlab 2) 
Kv = Cd*U*H0;                   % Vertical Eddy Diffusivity; this should be 
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
deltaT = 200;            % Time step [s]. Must satisfy Courant condition.
T = (12*60+25)*60;       % M2 and M4 tide. Time in seconds (same at Td2)
Tend = 10*T;             % Five tidal periods modeled -> for very fine
                         % sand and large erosion constants more tidal
                         % periods need to be solved
t = 0:deltaT:Tend;
Nt = length(t);

ampD1 = 0;               % D1 does not need to be prescribed here.
ampM2 = 1.5481;
ampM4 = 0.2;
phaseD1 = 0;
phaseM2 = 0.4352;        % Phase from Matlab 2
phaseM4 = 2.0975;        % Phase from Matlab 2

% Check if parameterisation obeys the Courant criterion
courant = sqrt(9.8*max(H0))*deltaT/deltaX;

if courant < 1
    D1 = ['Courant no. = ', num2str(courant)];
    disp(D1);
    disp('Courant criterion has been met');
end 

if courant >= 1
    D2 = ['Courant no. = ', num2str(courant)];
    disp(D2);
    disp('Courant criterion not met');
end 

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

% Assign depth based on equilibrium depth                         
H = H0*ones(1,length(x));

% Change in depth
dHdx(1:Nx) = 0;
% Alternative depth formulation for a linearly-sloping bottom.
% H = 10 + dHdx(1)*x; 

%% Find the flow velocity for all t, x
U = HydroModel2(t, Z, dZdt, H, dHdx, x, dx);

%% Find the sediment concentration for all t, x
for px = 1:Nx
    [C(px, 1:Nt)] = GroenModel(U(px, 1:Nt), t, deltaT, T, Ws, alpha, Kv);
end

%% Find the mean sediment transport for all x

% First find the sediment flux Qs ...
Qs = U.*C;

% ... then calculate the tidally-averaged sediment transport (averaged
% over the last tidal cycle)
S_Qs = 0;
for time = 1788:2012
    S_Qs = S_Qs + Qs(time);
end
meanQs = S_Qs/224; 

% Calculate tidally averaged sediment transport as a function of position
% in the estuary (only averaging over last tidal cycle)

Qs_x = [];
U_x = [];
for position = 1:81
    Qs_t(position)=0;
    U_t(position)=0;
for time = 1788:2012
    Qs_t(position) = Qs_t(position) + Qs(position,time);
    U_t(position) = U_t(position) + U(position,time);
end
% Tidally-averaged sediment transport
meanQs_x = Qs_t(position)/224;   
Qs_x = [Qs_x meanQs_x];

meanU_x = U_t(position)/224;
U_x = [U_x meanU_x];
end

%% Harmonic Analysis
% Define frequencies to be analysed. M1 tide is not included.
% wn is used by the harmfit function.
global wn
wn(1)=2*pi/Td2;         % M2
wn(2)=2*wn(1);          % M4
wn(3)=3*wn(1);          % M6

%Harmonic analysis
for px = 1:Nx
    coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
    coefout = nlinfit(t, U(px,:), @harmfit, coefin);
    U0(px) = coefout(1);
    UM2(px) = sqrt(coefout(2).^2+coefout(5).^2);
    UM4(px) = sqrt(coefout(3).^2+coefout(6).^2);
    UM6(px) = sqrt(coefout(4).^2+coefout(7).^2);
    phaseUM2(px) = atan(coefout(2)/coefout(5));
    phaseUM4(px) = atan(coefout(3)/coefout(6));
    phaseUM6(px) = atan(coefout(4)/coefout(7));
end
%% Plots 
phaseDiff = (2*phaseUM2(10) - phaseUM4(10))/pi*180

%Plot relation of tidally averaged sediment transport to spatial pattern of mean flow 
figure
subplot(2,1,1);
yyaxis right
plot(t(1788:2012),U(20,1788:2012))
ylabel('Flow velocity [m/s]');
xlabel('Time [s]');
hold on
plot(t(1788:2012),U(40,1788:2012))
plot(t(1788:2012),U(60,1788:2012))
plot(t(1788:2012),U(80,1788:2012))
hold off
yyaxis left
plot(t(1788:2012),Qs(20,1788:2012))
ylabel('Flux [kg m^{-1} s^{-1}]');
hold on
plot(t(1788:2012),Qs(40,1788:2012))
plot(t(1788:2012),Qs(60,1788:2012))
plot(t(1788:2012),Qs(80,1788:2012))
hold off
legend('40km','80km','120km','160km');
title('Sediment transport and flow velocity over time for multiple positions in the esturary');

subplot(2,1,2);
yyaxis left
plot(x/1000,Qs_x)
ylabel('Flux [kg m^{-1} s^{-1}]');
hold on 
yyaxis right
plot(x/1000,U_x)
plot(x/1000,UM2)
plot(x/1000,UM4)
ylabel('Flow velocity [m/s]');
xlabel('x [km]');
legend('Tidally averaged sediment transport Qs_U','TIdally averaged flow velocity U','Tidally averaged M2 flow U_{M2}','Tidally averaged M4 flow U_{M4}');
title('Tidally averaged flow velocities and sediment transport')
grid on;
hold off
savefig('Matlab3_3_i');

figure
subplot(2,1,1);
plot(U_x,Qs_x);
xlabel('Velocity [m/s]');
ylabel('Flux [kg m^{-1} s^{-1}]');
title('Sediment transport vs. average flow');
grid on,

subplot(2,1,2);
plot(UM2, Qs_x);
hold on
plot(UM4, Qs_x);
xlabel('Velocity [m/s]');
hold off
ylabel('Flux [kg m^{-1} s^{-1}]');
title('Sediment transport vs. velocity of tidal component');
legend('U_{M2}','U_{M4}');
grid on;
savefig('Matlab3_3_ii');
%% Mean sediment transport due to ... NEW FORMULA 
% presence of mean flows
qb_mf = alpha*U0.^3 + (3*alpha/2).*(UM2.^2).*U0;

% tidal asymmetry
qb_ta = (3*alpha/4).*(UM2.^2).*UM4.*cos(2*phaseM2 - phaseM4);

% Plot of both of the above
figure
plot(x/1000,qb_mf*1000);
hold on
plot(x/1000,qb_ta*1000);
hold off
xlabel('x [km]');
ylabel('qs [mms^{-1}]');
grid on;
legend('mean flow','tidal asymmetry');
title('Sediment transport by mean flow and tidal asymmetry');
savefig('Matlab3_3_iv');

%% Answers to Part 3.
% This will need to be rewritten/updated.
% Figure 1. The mean flow appears to entrain significant amounts of sediment.
% Tidal asymmetry accounts for the net difference in sediment transport.
% More specifically, we can say that duration asymmetry is not relevant for
% the sediment transport. Velocity asymmetry, however, is critical. We can 
% see that since the flood velocity is larger, the flood current entrains
% greater amounts of sediment, which means there will be a net flow of 
% sediment in the flood direction. 
% Figure 2. We can see a net flow in the flood direction on the order of 
% mm/s. Velocities and sediment concentration decline as we move towards
% the landward of the basin (as we would expect). 