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
global wn
wn(1)=2*pi/T;
wn(2)=2*wn(1);

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
%PhaseM4=pi/2; 

%sensitivity analysis phase difference
PhaseM4=(0:90:90)/180*pi;
Np=length(PhaseM4);

DIFF_phase=[];
figure
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
%Plot to analyse sensitivity of sediment concentration and transport to
%advection
if i==1
subplot(2,1,1);
title('Sediment Concentration and Transport when phaseM2-phaseM4=0');
elseif i==2  
subplot(2,1,2);
title('Sediment Concentration and Transport when phaseM2-phaseM4={\pi/2}');
end
yyaxis left
plot(t(235:418)/3600,C(1,235:418)*1000)
hold on
plot(t(235:418)/3600,C2(1,235:418)*1000)
ylabel('Concentration [g/m^2]');
yyaxis right
plot(t(235:418)/3600,Qs(1,235:418)*1000)
plot(t(235:418)/3600,Qs2(1,235:418)*1000)
ylabel('Flux [g/(m*s)]');
xlabel('Time [h]');
legend('Concentration without advection','Concentration with advection','Transport without advection','Transport with advection');
hold off
savefig('Matlab3_2_i');

%figure
%plot(t,U)

%Sediment concentration at flood and ebb (g/m^2)
C_Max(i,:)=findpeaks(C(1,285:418))*1000;
C2_Max(i,:)=findpeaks(C2(1,285:418))*1000;

%Sediment transport at flood and ebb (g/m*s)
Qs_Max(i)=findpeaks(Qs(1,235:335))*1000;
Qs2_Max(i)=findpeaks(Qs2(1,235:335))*1000;
Qs_Min(i)=findpeaks(-Qs(1,235:335))*-1000;
Qs2_Min(i)=findpeaks(-Qs2(1,235:335))*-1000;

% DIFF_S=[];
% for q=1:2
%     Diff=C_Max(i,q)-C2_Max(i,q);
%     DIFF_S=[DIFF_S Diff];
% end
% display(DIFF_S);

% %Difference between sediment concentration at flood and ebb (g/m^2)
DIFF_S2(i,:)=[C_Max(i,1)-C2_Max(i,1) C_Max(i,2)-C2_Max(i,2)]

%Difference between sediment transport at flood and ebb (g/m*s)
DIFF_Qs(i,:)=[Qs_Max(i)-Qs2_Max(i) Qs_Min(i)-Qs2_Min(i)]
end
%% 2.1.b
clear all
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
global wn
wn(1)=2*pi/T;
wn(2)=2*wn(1);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
%Remove amplitude of M4
ampM4=0;
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
for px=1:Nx
coefin=[0.1, 0.1, 0.1, 0.1, 0.1];
coefout=nlinfit(t,U(px,:),@harmfit,coefin);
U0(px)=coefout(1);
UM2(px)=sqrt(coefout(2).^2+coefout(3).^2);
UM4(px)=sqrt(coefout(3).^2+coefout(5).^2);
phaseUM2(px)=atan(coefout(2)/coefout(3));
phaseUM4(px)=atan(coefout(3)/coefout(5));
end

UM2=UM2';
Qs_M2=UM2.*C;                                            % Qs is sediment flux

% calculate tidally averaged sediment transport for M2 signal (only averaging over last tidal cycle)
S_Qs_M2=0;
for time=1342:1491
    S_Qs_M2=S_Qs_M2+Qs_M2(time);
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
    [C2(px,1:Nt)]=C1Model(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv, dQsdx(px,1:Nt));
end

%Sediment transport with advection
Qs2=U.*C2;

% calculate tidally averaged sediment transport with advection (only averaging over last tidal cycle)
S_Qs2=0;
for time=1342:1491
    S_Qs2=S_Qs2+Qs2(time);
end

% tidally averaged sediment transport
meanQs2=S_Qs2/149;

%Sediment transport with advection without M4
Qs2_M2=UM2.*C2;  

% calculate tidally averaged sediment transport for M2 signal (only averaging over last tidal cycle)
S_Qs2_M2=0;
for time=1342:1491
    S_Qs2_M2=S_Qs2_M2+Qs2_M2(time);
end

% tidally averaged sediment transport
meanQs2_M2=S_Qs2_M2/149; 

% Calculate the variation in tidally averaged sediment transport with and
% without M4 signal for the model with advection and the model without
% advection. 

DIFF_Qs_M4=meanQs-meanQs_M2;
DIFF_Qs2_M4=meanQs2-meanQs2_M2;

%Plot of difference in transport and flux as a function of space TO BE DONE
% figure
% ylabel('Flux difference [kg/(m*s)]');
% plot(x,abs(Qs2(:,1342:1491)-Qs2_M2))
% hold on
% scatter(abs(DIFF_QS2_M4))
% xlabel('Ws [m/s]');
% legend('Difference in sediment transport without advection','Difference in sediment transport with advection');
% title('Sensiivity Analysis of sediment transport for varying fall velocities (Ws [m/s])')
% hold off
% savefig('Matlab3_2_ii');

%% 2.1.c 
clear all
%Ws=1e-3;                % Fall velocity of sediment
alpha=1e-4;             % Erosion coefficent
Kv=1e-2;                % Vertical eddy diffusivity (for vertical mixing)

%Sensitivity analysis Ws
WS=linspace(0.5e-3,2e-2,5);

DIFF_QS_M4 = [];
DIFF_QS2_M4 = [];
for j=1:5
    Ws=WS(j);

%**************************************************************************
%           Define time domain
%**************************************************************************

T=(12*60+25)*60;        % We model only the M2 and M4 tide. Time is in seconds. 
Tend=10*T;               % Five tidal periods modeled -> for very fine sand and large erosion constants more tidal periods need to be solved
deltaT=300;             % Time step of 5 minutes
t=0:deltaT:Tend;
Nt=length(t);
global wn
wn(1)=2*pi/T;
wn(2)=2*wn(1);

%**************************************************************************
% Prescribed sea surface elevations. It is assumed that d/dx zeta =0. M2
% and M4 are prescribed at the seaward boundary. So these are the sea surface heights in the entire
% basin at any moment in time.
%**************************************************************************
ampD1=0;            % in part 1 and 2 D1=0. Depending on your estuary, you might want to prescribe D1 for part 3. 
ampM2=1;
%Remove amplitude of M4
ampM4=0;
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
for px=1:Nx
coefin=[0.1, 0.1, 0.1, 0.1, 0.1];
coefout=nlinfit(t,U(px,:),@harmfit,coefin);
U0(px)=coefout(1);
UM2(px)=sqrt(coefout(2).^2+coefout(3).^2);
UM4(px)=sqrt(coefout(3).^2+coefout(5).^2);
phaseUM2(px)=atan(coefout(2)/coefout(3));
phaseUM4(px)=atan(coefout(3)/coefout(5));
end

if j==1
UM2=UM2';
end
Qs_M2=UM2.*C;                                            % Qs is sediment flux

% calculate tidally averaged sediment transport for M2 signal (only averaging over last tidal cycle)
S_Qs_M2=0;
for time=1342:1491
    S_Qs_M2=S_Qs_M2+Qs_M2(time);
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
    [C2(px,1:Nt)]=C1Model(U(px,1:Nt),t,deltaT, T, Ws, alpha, Kv, dQsdx(px,1:Nt));
end

%Sediment transport with advection
Qs2=U.*C2;

% calculate tidally averaged sediment transport with advection (only averaging over last tidal cycle)
S_Qs2=0;
for time=1342:1491
    S_Qs2=S_Qs2+Qs2(time);
end

% tidally averaged sediment transport
meanQs2=S_Qs2/149;

%Sediment transport with advection without M4
Qs2_M2=UM2.*C2;  

% calculate tidally averaged sediment transport for M2 signal (only averaging over last tidal cycle)
S_Qs2_M2=0;
for time=1342:1491
    S_Qs2_M2=S_Qs2_M2+Qs2_M2(time);
end

% tidally averaged sediment transport
meanQs2_M2=S_Qs2_M2/149; 

% Calculate the variation in tidally averaged sediment transport with and
% without M4 signal for the model with advection and the model without
% advection. 

DIFF_Qs_M4=meanQs-meanQs_M2;
DIFF_Qs2_M4=meanQs2-meanQs2_M2;

% Save DIFF_Qs_M4 and DIFF_Qs2_M4 for each value of WS.    
DIFF_QS_M4(j) = DIFF_Qs_M4;
DIFF_QS2_M4(j) = DIFF_Qs2_M4;
end

%Plot of sensitivity analysis
figure
plot(WS,abs(DIFF_QS2_M4)*1000)
ylabel('Flux difference [g/(m*s)]');
hold on
scatter(WS,abs(DIFF_QS2_M4)*1000)
xlabel('Ws [m/s]');
%legend('Difference in sediment transport without advection','Difference in sediment transport with advection');
title('Sensiivity Analysis of sediment transport with advection for varying fall velocities (Ws [m/s])')
hold off
savefig('Matlab3_2_ii');


