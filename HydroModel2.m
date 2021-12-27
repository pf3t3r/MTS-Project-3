function U=HydroModel2(t,Z,dZdt,H,dHdx,x,dx)


% Hydromodel.
% t = time
% Z is water level at seaward boundary
% dZdt change of water level in time
% H is still water depth as a function of position x
% L is length of basin
% dx is griddistance


Nx=length(x);
Nt=length(t);
L=x(end);

testPlane=cumsum(abs(dHdx)>0);

% solve dZdt + d/dx Hu=0, assuming dZdx=0. So U*dHdx + (H+Z)*dU/dx = -dZ/dt
% U=0 at seaward boundary.

if testPlane==0
    H0=H(1);
    U(1:Nx,1:Nt)=x'*(-dZdt./(H0+Z));
else
    U(Nx,1:Nt)=0.0;
    for px=1:Nx-1
        U(px,1:Nt)=(x(px)-L)*(-dZdt./(H(px)+Z));
    end
end

