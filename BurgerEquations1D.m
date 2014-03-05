% AA543 - CP 3.4
% Tyler James Pierce
% March 04 2014
% Simulating the 1-D Burger equation

problemNumber = 1; %Choose which problem to run

% Parameters
nx=30;                          %number of steps in space(x)
nt=60;                           %number of time steps 
dt=0.01;                         %Width of each time step
rangex=[-pi,pi];
L=max(rangex)-min(rangex);%
vis=0.01;                        %Diffusion coefficient/viscosity
u=zeros(1,nx);                  %Preallocating u
un=zeros(1,nx);                 %Preallocating un
rhsx=zeros(1,nx);
dxuniform=L/(nx-1);

switch problemNumber
  case 1
    a = 1; %amplitude
    BC = 'periodic';
    x = rangex(1):dxuniform:rangex(end);
  case 2
    a = -1;     %amplitude
    BC = 'dirichlet';
    x = rangex(1):dxuniform:rangex(end);
  case 3
    a = -1; %amplitude
    SF = 0.8; %grid clustering
    BC = 'dirichlet';
    xc = 0;
    x = Mesh1D(range,nx,xc,SF);
end
u0=a*sin(x);
u=u0;

