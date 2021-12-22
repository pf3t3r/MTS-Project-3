clear
clc
close all

%**************************************************************************
% Parameters needed for the sediment transport calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 1.a.i
%Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

% Sensitivity analysis Ws
WS=linspace(0.5e-3,2e-2,5);

DIFF_ws=[];
for i=1:5
    Ws=WS(i);
%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
ampM4=0.2;
phaseD1=0;
phaseM2=0;
phaseM4=pi/2;

Z=ampD1*sin(pi*t/T + phaseD1)+ampM2*sin(2*pi*t/T + phaseM2)+ampM4*sin(4*pi*t/T + phaseM4);          % Waterlevel prescribed as sine function. 
dZdt=ampD1*1*pi/T*cos(pi*t/T+ phaseD1)+ampM2*2*pi/T*cos(2*pi*t/T+ phaseM2)+ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); % Flow velocity will behave as a cosine function. 

%**************************************************************************
%       Spatial Domain and Grid
%**************************************************************************

L=1e4;                  % We model a simple basin with a length of ten km
dx=400;                 % Grid distance
x=0:dx:L;               % x-coordinate. Seaward end is at x=L, landward end at x=0. 
Nx=length(x);                       
                        
%**************************************************************************
%
%   x=0 (=Inlet) ...................... x=L (=landward side of basin)
%
%   So x=positive in landward direction
%
%   U>0 = Flood flow          U<0 = Ebb flow
%
%**************************************************************************

%**************************************************************************
%           Bed level in basin
%**************************************************************************

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom. 2 m deep near landward boundary, 10 m deep near inlet. 
dHdx(1:Nx)=-8e-4;

%**************************************************************************
% After a call to hydromodel flow velocity at each position as a function of
% time is known
%**************************************************************************

U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

%**************************************************************************
% Here you have to calculate the sediment concentrations with the Groen
% model. This is a Matlab function which has as input the flow velocity, the relevant
% parameters, and time. For each position in the basin do a call to this
% Groenmodel. 

for px=1:Nx
    [C(px,1:Nt)]=GroenModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv);
end

%Spot the correct section to consider for U_Max and C_Max
% figure
% yyaxis left
% plot(t,C(1,:))
% hold on
% yyaxis right
% plot(t,U(1,:))
% hold off

%Peak sediment concentration and time 
C_Max=max(C(1,101:201));
t_C_Max=t(C(1,:)==C_Max);

%Peak flood flow and time 
U_Max=max(U(1,101:201));
t_U_Max=t(U(1,:)==U_Max);

%Time difference between peaks
Diff=(t_C_Max-t_U_Max)/3600;
DIFF_ws=[DIFF_ws Diff];
end 

figure
plot(WS,abs(DIFF_ws))
title('Sensitivity analysis of time difference between peak flow and peak sediment concentration for varying fall velocities (W_s)');
xlabel('W_{s} [m/s]');
ylabel('Concentration [kg/m^2]');
grid on;
savefig('Matlab3_1_i');

%% Part 1.a.ii

Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
%Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

% Sensitivity analysis Kv
KV=linspace(1e-3,1e-1,5);

DIFF_kv=[];
for i=1:5
    Kv=KV(i);
%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
ampM4=0.2;
phaseD1=0;
phaseM2=0;
phaseM4=pi/2;

Z=ampD1*sin(pi*t/T + phaseD1)+ampM2*sin(2*pi*t/T + phaseM2)+ampM4*sin(4*pi*t/T + phaseM4);          % Waterlevel prescribed as sine function. 
dZdt=ampD1*1*pi/T*cos(pi*t/T+ phaseD1)+ampM2*2*pi/T*cos(2*pi*t/T+ phaseM2)+ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); % Flow velocity will behave as a cosine function. 

%**************************************************************************
%       Spatial Domain and Grid
%**************************************************************************

L=1e4;                  % We model a simple basin with a length of ten km
dx=400;                 % Grid distance
x=0:dx:L;               % x-coordinate. Seaward end is at x=L, landward end at x=0. 
Nx=length(x);                       
                        
%**************************************************************************
%
%   x=0 (=Inlet) ...................... x=L (=landward side of basin)
%
%   So x=positive in landward direction
%
%   U>0 = Flood flow          U<0 = Ebb flow
%
%**************************************************************************

%**************************************************************************
%           Bed level in basin
%**************************************************************************

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom. 2 m deep near landward boundary, 10 m deep near inlet. 
dHdx(1:Nx)=-8e-4;

%**************************************************************************
% After a call to hydromodel flow velocity at each position as a function of
% time is known
%**************************************************************************

U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

%**************************************************************************
% Here you have to calculate the sediment concentrations with the Groen
% model. This is a Matlab function which has as input the flow velocity, the relevant
% parameters, and time. For each position in the basin do a call to this
% Groenmodel. 

for px=1:Nx
    [C(px,1:Nt)]=GroenModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv);
end


%Peak sediment concentration and time 
C_Max=max(C(1,101:201));
t_C_Max=t(C(1,:)==C_Max);

%Peak flood flow and time 
U_Max=max(U(1,101:201));
t_U_Max=t(U(1,:)==U_Max); 

%Time difference between peaks
Diff=(t_C_Max-t_U_Max)/3600;
DIFF_kv=[DIFF_kv Diff];
end 

figure
plot(KV,abs(DIFF_kv))
title('Sensitivity analysis of time difference between peak flow and peak sediment concentration for varying eddy diffusivities (K_v)');
xlabel('K_{v} [m^{2}/s]');
ylabel('Concentration [kg/m^2]');
grid on;
savefig('Matlab3_1_ii');

%% Part 1.b.i
Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
ampM4=0.2;
phaseD1=0;
phaseM2=0;
%phaseM4=pi/2;

%sensitivity analysis phase difference
PhaseM4=(0:45:180)/180*pi;
Np=length(PhaseM4);

DIFF_phase=[];
for i=1:Np
    phaseM4=PhaseM4(i);
Z=ampD1*sin(pi*t/T + phaseD1)+ampM2*sin(2*pi*t/T + phaseM2)+ampM4*sin(4*pi*t/T + phaseM4);          % Waterlevel prescribed as sine function. 
dZdt=ampD1*1*pi/T*cos(pi*t/T+ phaseD1)+ampM2*2*pi/T*cos(2*pi*t/T+ phaseM2)+ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); % Flow velocity will behave as a cosine function. 

%**************************************************************************
%       Spatial Domain and Grid
%**************************************************************************

L=1e4;                  % We model a simple basin with a length of ten km
dx=400;                 % Grid distance
x=0:dx:L;               % x-coordinate. Seaward end is at x=L, landward end at x=0. 
Nx=length(x);                       
                        
%**************************************************************************
%
%   x=0 (=Inlet) ...................... x=L (=landward side of basin)
%
%   So x=positive in landward direction
%
%   U>0 = Flood flow          U<0 = Ebb flow
%
%**************************************************************************

%**************************************************************************
%           Bed level in basin
%**************************************************************************

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom. 2 m deep near landward boundary, 10 m deep near inlet. 
dHdx(1:Nx)=-8e-4;

%**************************************************************************
% After a call to hydromodel flow velocity at each position as a function of
% time is known
%**************************************************************************

U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

%**************************************************************************
% Here you have to calculate the sediment concentrations with the Groen
% model. This is a Matlab function which has as input the flow velocity, the relevant
% parameters, and time. For each position in the basin do a call to this
% Groenmodel. 

for px=1:Nx
    [C(px,1:Nt)]=GroenModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv);
end

%Tracking the behaviours and locating the peaks
% figure
% yyaxis left
% plot(t,C(1,:))
% hold on
% yyaxis right
% plot(t,U(1,:))
% hold off

%Peak sediment concentration at flood and ebb
C_Max=findpeaks(C(1,201:371));

%Difference between peaks
if length(C_Max)==1
    Diff=abs(C_Max(1)-max(C(1,268:371)));
else
    Diff=C_Max(1)-C_Max(2);
end
DIFF_phase=[DIFF_phase Diff];
end

figure
plot(PhaseM4,abs(DIFF_phase))
title('Sensitivity analysis of difference in peak sediment concentration between peak flow and peak ebb for varying relative phase differences between M2 and M4');
xlabel('Phase difference [rad]');
ylabel('Concentration [kg/m^2]');
grid on;
savefig('Matlab3_1_iii');

%% Part 1.b.ii
%Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

% Sensitivity analysis Ws
WS=linspace(0.5e-3,2e-2,5);
 
DIFF_ws_ii=[];
for i=1:5
    Ws=WS(i);

%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
ampM4=0.2;
phaseD1=0;
phaseM2=0;
phaseM4=pi/2;

Z=ampD1*sin(pi*t/T + phaseD1)+ampM2*sin(2*pi*t/T + phaseM2)+ampM4*sin(4*pi*t/T + phaseM4);          % Waterlevel prescribed as sine function. 
dZdt=ampD1*1*pi/T*cos(pi*t/T+ phaseD1)+ampM2*2*pi/T*cos(2*pi*t/T+ phaseM2)+ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); % Flow velocity will behave as a cosine function. 

%**************************************************************************
%       Spatial Domain and Grid
%**************************************************************************

L=1e4;                  % We model a simple basin with a length of ten km
dx=400;                 % Grid distance
x=0:dx:L;               % x-coordinate. Seaward end is at x=L, landward end at x=0. 
Nx=length(x);                       
                        
%**************************************************************************
%
%   x=0 (=Inlet) ...................... x=L (=landward side of basin)
%
%   So x=positive in landward direction
%
%   U>0 = Flood flow          U<0 = Ebb flow
%
%**************************************************************************

%**************************************************************************
%           Bed level in basin
%**************************************************************************

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom. 2 m deep near landward boundary, 10 m deep near inlet. 
dHdx(1:Nx)=-8e-4;

%**************************************************************************
% After a call to hydromodel flow velocity at each position as a function of
% time is known
%**************************************************************************

U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

%**************************************************************************
% Here you have to calculate the sediment concentrations with the Groen
% model. This is a Matlab function which has as input the flow velocity, the relevant
% parameters, and time. For each position in the basin do a call to this
% Groenmodel. 

for px=1:Nx
    [C(px,1:Nt)]=GroenModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv);
end

%Tracking the behaviours and locating the peaks
% figure
% yyaxis left
% plot(t,C(1,:))
% hold on
% yyaxis right
% plot(t,U(1,:))
% hold off

%Peak sediment concentration at flood and ebb
C_Max=findpeaks(C(1,201:371));

%Difference between peaks
if length(C_Max)==1
    Diff=abs(C_Max(1)-max(C(1,268:371)));
else
    Diff=C_Max(1)-C_Max(2);
end
DIFF_ws_ii=[DIFF_ws_ii Diff];
end

figure
plot(WS,abs(DIFF_ws_ii))
title('Sensitivity analysis of difference in peak sediment concentration between peak flow and peak ebb for varying fall velocities (W_s)');
xlabel('W_{s} [m/s]');
ylabel('Concentration [kg/m^2]');
grid on;
savefig('Matlab3_1_iv');

%% Part 1.b.iii
Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
%Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

% Sensitivity analysis Kv
KV=linspace(1e-3,1e-1,5);

DIFF_kv_ii=[];
for i=1:5
    Kv=KV(i);

%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
ampM4=0.2;
phaseD1=0;
phaseM2=0;
phaseM4=pi/2;

Z=ampD1*sin(pi*t/T + phaseD1)+ampM2*sin(2*pi*t/T + phaseM2)+ampM4*sin(4*pi*t/T + phaseM4);          % Waterlevel prescribed as sine function. 
dZdt=ampD1*1*pi/T*cos(pi*t/T+ phaseD1)+ampM2*2*pi/T*cos(2*pi*t/T+ phaseM2)+ampM4*4*pi/T*cos(4*pi*t/T + phaseM4); % Flow velocity will behave as a cosine function. 

%**************************************************************************
%       Spatial Domain and Grid
%**************************************************************************

L=1e4;                  % We model a simple basin with a length of ten km
dx=400;                 % Grid distance
x=0:dx:L;               % x-coordinate. Seaward end is at x=L, landward end at x=0. 
Nx=length(x);                       
                        
%**************************************************************************
%
%   x=0 (=Inlet) ...................... x=L (=landward side of basin)
%
%   So x=positive in landward direction
%
%   U>0 = Flood flow          U<0 = Ebb flow
%
%**************************************************************************

%**************************************************************************
%           Bed level in basin
%**************************************************************************

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom. 2 m deep near landward boundary, 10 m deep near inlet. 
dHdx(1:Nx)=-8e-4;

%**************************************************************************
% After a call to hydromodel flow velocity at each position as a function of
% time is known
%**************************************************************************

U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx);

%**************************************************************************
% Here you have to calculate the sediment concentrations with the Groen
% model. This is a Matlab function which has as input the flow velocity, the relevant
% parameters, and time. For each position in the basin do a call to this
% Groenmodel. 

for px=1:Nx
    [C(px,1:Nt)]=GroenModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv);
end

%Tracking the behaviours and locating the peaks
% figure
% yyaxis left
% plot(t,C(1,:))
% hold on
% yyaxis right
% plot(t,U(1,:))
% hold off

%Peak sediment concentration at flood and ebb
C_Max=findpeaks(C(1,201:371));

%Difference between peaks
if length(C_Max)==1
    Diff=abs(C_Max(1)-max(C(1,268:371)));
else
    Diff=C_Max(1)-C_Max(2);
end
DIFF_kv_ii=[DIFF_kv_ii Diff];
end

figure
plot(kV,abs(DIFF_kv_ii))
title('Sensitivity analysis of difference in peak sediment concentration between peak flow and peak ebb for varying eddy diffusivities (K_v)');
xlabel('K_{v} [m^{2}/s]');
ylabel('Concentration [kg/m^2]');
grid on;
savefig('Matlab3_1_v');

%Qs=U.*C;                                            % Qs is sediment flux

%Nsteps=T/deltaT;                                    % Nr of timestepf in one tidal cycle.

%Nsteps=2*T/deltaT                                  % If you have precribed a D1 tide

% calculate tidally averaged sediment transport (only averaging over last tidal cycle)
%
%           meanQs=
%           
%**************************************************


%**************************************************************************
% Use calculated flux with Groen's model to calculate sediment
% concentration with model that includes advective fluxes
%
% Make Cmodel yourself.
% It should include a d/dx UC term
% You can use code like this
% [dQsdt dQsdx]=gradient(Qs,deltaT,dx);
%**************************************************************************

% for px=1:Nx
%     [C2(px,1:Nt)]=CModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv, dQsdx(px,1:Nt));
% end
% 
% Qs2=U.*C2;


