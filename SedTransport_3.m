clc; clear; close all;

%% Analysis of Sediment Transport in the Gironde Estuary

% Coefficients and typical values
alpha = 1e-4;                   % Erosion coefficent
Ws = 3e-3;                      % Fall velocity
Cd = 0.2e-3;                    % Typical drag coefficient [] (Matlab 2)
U = 8;                          % Typical flow velocity [m/s] (Matlab 2)
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
% time = 0:deltaT:30*Td2;  % Time [s] --> already defined below!!
% Nt = length(time);       % Size of time array

T = (12*60+25)*60;        % M2 and M4 tide. Time in seconds (same at Td2)
Tend = 10*T;              % Five tidal periods modeled -> for very fine
                          % sand and large erosion constants more tidal
                          % periods need to be solved
% deltaT = 300;             % Time step of 5 minutes -- already defined!!!
t = 0:deltaT:Tend;
Nt = length(t);

ampD1 = 0;                % D1 does not need to be prescribed here.
ampM2 = 1;
ampM4 = 0.2;
phaseD1 = 0;
phaseM2 = 0;
phaseM4 = 60/180*pi;      % Phase = 60 deg. Might not be correct

% Check if parameterisation obeys Courant 
courant = sqrt(9.8*max(H0))*deltaT/deltaX;

if courant<1
    D1 = ['Courant no. = ', num2str(courant)];
    disp(D1);
    disp('Courant criterion has been met');
end 

if courant>=1
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

% Not sure if the following is correct.
% H = 10-8e-4*x;         % Bottom profile. Linear sloping bottom. 2 m deep
                         % near landward boundary, 10 m deep near inlet. 
H = H0*ones(1,length(x));
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

Qs=U.*C;                                            % Qs is sediment flux

Nsteps=T/deltaT;                                    % Nr of timestepf in one tidal cycle.

% calculate tidally averaged sediment transport (only averaging over last tidal cycle)
S_Qs=0;
for time=2013:2236
    S_Qs=S_Qs+Qs(time);
end

% tidally averaged sediment transport
meanQs=S_Qs/223; 

% calculate tidally averaged sediment transport as a function of position in the estuary (only averaging over last tidal cycle)
Qs_x=[];
for position=1:81
Qs_t(position)=0;
for time=2013:2236
    Qs_t(position)=Qs_t(position)+Qs(position,time);
end
%tidally averaged sediment transport
meanQs_x=Qs_t(position)/223;   
Qs_x=[Qs_x meanQs_x];
end
%**************************************************

%Plot of tidally averaged sediment transport as a funcion of position in
%the estuary 
figure
ylabel('Flux [kg/(m*s)]');
plot(x,Qs_x)
xlabel('x [m]');
title('Tidally averaged sediment transport as a function of position in the estuary')
savefig('Matlab3_3_i');
%% Determine tidally-averaged sediment transport as a function of position
% (and also tidally-averaged velocity)
for px = 1:Nx
    C_avg(px) = mean(C(px,:));
    U_avg(px) = mean(U(px,:));
end

% Comment the following to suppress ouput
disp('C_avg: ')
disp(C_avg)

% Create a legend
 for i = 1:Nx
      C_legend{i} = num2str(x(i)/1000,'x = %.0f km');
      U_legend{i} = num2str(x(i)/1000,'x = %.0f km');
 end

% Decide how many x locations will be shown on the time plot
% This is mostly just a workaround for the fact that the fourth line
% that is automatically drawn has circles and looks like shit.
xlocs = 3;

% Sediment Concentration & Velocity VS. Time
figure

subplot(2,1,1)
yyaxis left
for i = 1:int16(Nx/xlocs):Nx
    plot(t/3600,C(i,:));
    hold on
end
ylabel('C [kg/m^2]');
legend(C_legend(1:int16(Nx/xlocs):Nx));
yyaxis right
for i = 1:int16(Nx/xlocs):Nx
    plot(t/3600,U(i,:));
end
legend(U_legend(1:int16(Nx/xlocs):Nx));
hold off
ylabel('U [m/s]');
xlabel('t [hrs]');
grid(gca,'minor')
grid on;
title('U, C vs. t')

subplot(2,1,2)
yyaxis left
for i = 1:int16(Nx/xlocs):Nx
    plot(t(1:500)/3600,C(i,1:500));
    hold on
end
ylabel('C [kg/m^2]');
yyaxis right
for i = 1:int16(Nx/xlocs):Nx
    plot(t(1:500)/3600,U(i,1:500));
end
hold off
ylabel('U [m/s]');
xlabel('t [hrs]');
grid(gca,'minor')
grid on;
legend(U_legend(1:int16(Nx/xlocs):Nx));
title('U, C vs. t: first two tidal cycles');

savefig('pt-3-i');

% (Tidally-averaged) Sediment Concentration & Velocity VS. Basin Length

figure
yyaxis left
plot(x/1000,U_avg*1000);
ylabel('U [mm/s]');
hold on
yyaxis right
plot(x/1000,C_avg);
ylabel('C [kg/m^2]');
hold off
grid(gca,'minor')
grid on;
xlabel('Estuary Length [km]');
legend('U_{avg}','C_{avg}');
title('Tidally-averaged U, C vs. X');

savefig('pt-3-ii');

% The above shows the dependence of the tidally-averaged sediment transport
% on mean flow and tidal asymmetry. We still need to find the following.

%% Relationship of tidally-averaged sediment transport to magnitude of tidal flows
% this is the individual components

%% Relationship of tidally-averaged sediment transport to generation of higher harmonics

% Define frequencies to be analysed. M1 tide is not included.
global wn
wn(1)=2*pi/Td2;         % M2
wn(2)=2*wn(1);          % M4
wn(3)=3*wn(1);          % M6

% I am not even using wn yet. LOL.
% OK. Let's add that stuff below.

for px=1:Nx-1
    coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
    coefout=nlinfit(t(1:end),Z(1:end),@harmfit,coefin);
    Z0(px)=coefout(1);
    ZM2(px)=sqrt(coefout(2).^2+coefout(5).^2);
    ZM4(px)=sqrt(coefout(3).^2+coefout(6).^2);
    ZM6(px)=sqrt(coefout(4).^2+coefout(7).^2);
    phaseZM2(px)=atan(coefout(2)/coefout(5));
    phaseZM4(px)=atan(coefout(3)/coefout(6));
    phaseZM6(px)=atan(coefout(4)/coefout(7));
    % coefin=[0.1, 0.3, 1, 0.2, 0.1, 0.2, 1, 0.2, 0.1];
    coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
    coefout=nlinfit(t(1:end),U(px,1:end),@harmfit,coefin);
%     U0(px)=coefout(1);
%     UM2(px)=sqrt(coefout(2).^2+coefout(5).^2);
%     UM4(px)=sqrt(coefout(3).^2+coefout(6).^2);
%     UM6(px)=sqrt(coefout(4).^2+coefout(7).^2);
%     phaseUM2(px)=atan(coefout(2)/coefout(5));
%     phaseUM4(px)=atan(coefout(3)/coefout(6));
%     phaseUM6(px)=atan(coefout(4)/coefout(7));
end

%% Relationship of tidally-averaged sediment transport to the phase difference between M2 and M4 flow velocities


%% Answers to Part 3.
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