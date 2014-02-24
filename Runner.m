function [index,x,dx,u,time] = Runner(range,imax,xPeak,xc,SF,BC,IC,tfinal,niter,cfl)
%RUNNER Main function, runs all other internal functions

% Initialize variables
uold=zeros(1,imax); % row index = time, column index = location
unew=zeros(1,imax);
uexact=zeros(1,imax);
u=zeros(niter,imax);

% get Mesh
[index,x,dx]=Mesh1D(range,imax,xc,SF);
% Initial Condition
% change to uold=InitialCondition(IC)
if strcmpi(IC,'square')==1
    uold=SquareWave(x,0,1); %umin=0, umax=1
elseif strcmpi(IC,'step')==1
    uold=StepWave(x,xPeak,SF);
elseif strcmpi(IC,'gauss')==1
    uold=GaussianWave(range,x,xPeak);
end
uexact=uold;
u(1,:)=uold;
t=0;
iter=1;
dt=0;
time=zeros(1,tfinal);

  while t <= tfinal-dt && iter < niter  %Advance the solution to t+dt
    dxgrid=max(dx);
    dt=cfl*dxgrid;
    iter=iter+1;
    [unew]=GetNewProperties(cfl,uold);
    t=t+dt;
    time(iter)=t;
    [unew] = BoundaryCondition(BC, unew, [0 0]);
    uold=unew;
    u(iter,:)=unew;
  end
end

% Author: Tyler James Pierce
% Tylerpierce644@gmail.com
% Created: Feb 07 2014

% Update History:
