clear; clc; close all;

%**************************************************************************
% Parameters needed for the sediment transport calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 2: patial settling lag mechanism in short basin

%% 2.1:

%% 2.1.a 
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
PhaseM4=(0:90:90)/180*pi;
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

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom.
                         % 2 m deep near landward boundary, 10 m deep near inlet. 
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

Qs=U.*C;                                            % Qs is sediment flux

Nsteps=T/deltaT;                                    % Nr of timestepf in one tidal cycle.

% calculate tidally averaged sediment transport (only averaging over last tidal cycle)
S_Qs=0;
for time=1342:1491
    S_Qs=S_Qs+Qs(time);
end

% tidally averaged sediment transport
meanQs=S_Qs/149;           
%**************************************************


%**************************************************************************
% Use calculated flux with Groen's model to calculate sediment
% concentration with model that includes advective fluxes
%
% Make Cmodel yourself.
% It should include a d/dx UC term
% You can use code like this
[dQsdt dQsdx]=gradient(Qs,deltaT,dx);
%[dQsdt dQsdx]=gradient(meanQs,deltaT,dx);
%**************************************************************************

%Sediment concentration with advection
for px=1:Nx
    [C2(px,1:Nt)]=CModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv, dQsdx(px,1:Nt));
end
% 
Qs2=U.*C2;
% 
%For inspection
figure
yyaxis left
plot(t,C(1,:))
hold on
plot(t,C2(1,:))
yyaxis right
plot(t,Qs(1,:))
plot(t,Qs2(1,:))
hold off

%Sediment concentration at flood and ebb
C_Max(i,:)=findpeaks(C(1,268:435));
C2_Max(i,:)=findpeaks(C2(1,268:435));

%Sediment transport at flood and ebb
Qs_Max(i)=findpeaks(Qs(1,235:335));
Qs2_Max(i)=findpeaks(Qs2(1,235:335));
Qs_Min(i)=findpeaks(-Qs(1,235:335))*-1;
Qs2_Min(i)=findpeaks(-Qs2(1,235:335))*-1;

% DIFF_S=[];
% %Difference between sediment concentration at flood and ebb
% for q=1:2
%     Diff=C_Max(i,q)-C2_Max(i,q);
%     DIFF_S=[DIFF_S Diff];
% end
% display(DIFF_S);
DIFF_S2(i,:)=[C_Max(i,1)-C2_Max(i,1) C_Max(i,2)-C2_Max(i,2)]

%Difference between sediment concentration at flood and ebb
DIFF_Qs(i,:)=[Qs_Max(i)-Qs2_Max(i) Qs_Min(i)-Qs2_Min(i)]
end

%% 2.1.b 
Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

% Sensitivity analysis Ws
%WS=linspace(0.5e-3,2e-2,5);
%for i=1:5
    %Ws=WS(i);

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
phaseM4=0;
    
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

H=10-8e-4*x;             % Bottom profile. Linear sloping bottom.
                         % 2 m deep near landward boundary, 10 m deep near inlet. 
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

%Sediment transport and tidally averaged sediment transport
Qs=U.*C;                                            % Qs is sediment flux

Nsteps=T/deltaT;                                    % Nr of timestepf in one tidal cycle.

% calculate tidally averaged sediment transport (only averaging over last tidal cycle)
S_Qs=0;
for time=1342:1491
    S_Qs=S_Qs+Qs(time);
end

% tidally averaged sediment transport
meanQs=S_Qs/149; 

Nsteps=floor(T/deltaT);
%Harmonic analysis
for px=1:Nx-1
coefin=[0.1, 0.1, 0.1, 0.1, 0.1];
coefout=nlinfit(t(end-Nsteps:end),Z(px,end-Nsteps:end),@harmfit,coefin);
Z0(i,px)=coefout(1);
ZM2(i,px)=sqrt(coefout(2).^2+coefout(4).^2);
ZM4(i,px)=sqrt(coefout(3).^2+coefout(5).^2);
phaseZM2(i,px)=atan(coefout(2)/coefout(4));
phaseZM4(i,px)=atan(coefout(3)/coefout(5));
coefin=[0.1, 0.1, 0.1, 0.1, 0.1];
coefout=nlinfit(t(end-Nsteps:end),U(px,end-Nsteps:end),@harmfit,coefin);
U0(px)=coefout(1);
UM2(i,px)=sqrt(coefout(2).^2+coefout(4).^2);
UM4(i,px)=sqrt(coefout(3).^2+coefout(5).^2);
phaseUM2(i,px)=atan(coefout(2)/coefout(4));
phaseUM4(i,px)=atan(coefout(3)/coefout(5));
end

Qs_M2=UM2.*C;                                            % Qs is sediment flux

% calculate tidally averaged sediment transport (only averaging over last tidal cycle)
S_Qs_M2=0;
for time=1342:1491
    S_Qs_M2=S_Qs_M2+Qs-M2(time);
end

% tidally averaged sediment transport
meanQs_M2=S_Qs_M2/149;           
%**************************************************


%**************************************************************************
% Use calculated flux with Groen's model to calculate sediment
% concentration with model that includes advective fluxes
%
% Make Cmodel yourself.
% It should include a d/dx UC term
% You can use code like this
[dQsdt dQsdx]=gradient(Qs,deltaT,dx);
%[dQsdt dQsdx]=gradient(meanQs,deltaT,dx);
%**************************************************************************

%Sediment concentration with advection
for px=1:Nx
    [C2(px,1:Nt)]=CModel(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv, dQsdx(px,1:Nt));
end
% 
Qs2=U.*C2;
% 
