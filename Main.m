%Tyler Pierce 
%AA543
%Computer Project 3.3
%% (B) Square Wave with Dirichlet BC
%IC: Square Wave
%BC: Dirichlet (u1)=u(end)=0
%CFL=0.8 
%tmax @ square wave exits domain
%dx=constant
% xc=N/A; BC=du/dx=0 at x=min,xmax
% algorithm = Lax-Wendroff
SF=1;
range=[0,pi];
imax=401;
xpeak=0;
xc=NaN;
BC='Dirichlet';
IC='Square';
tfinal=imax;
nmax=1000;
cfl=0.8;

[index,x,dx,u,t]=Runner(range,imax,xpeak,xc,SF,BC,IC,tfinal,nmax,cfl);

%Create figure
figure1 = figure('Name','figure1u0vi','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
xlim(axes1,[0 max(x)]);
ylim(axes1,[0 1.4]);
box(axes1,'on');
hold(axes1,'on');
xlabel(axes1,'x');
ylabel(axes1,'u(x)');
title(axes1,{'Square Function Initial Wave u0(x)'});
plot(x,u(1,:),x,u(100,:),x,u(200,:),x,u(300,:),x,u(400,:),x,u(450,:));
% figure2 = figure('Name','figure1u0vx','Color',[1 1 1]);
% axes2 = axes('Parent',figure2);
% xlim(axes2,[min(range) max(range)]);
% ylim(axes2,[min(u) max(u)]);
% box(axes2,'on');
% hold(axes2,'on');
% xlabel(axes2,{'x_i'});
% ylabel(axes2,{'u(x_i)'});
% title(axes2,{'Step Function Initial Wave u0(x_i)'});
% plot(x,u,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
%     'MarkerSize',10,'Marker','.','Color',[0 0 0]);

%% (E) Square Wave @ CFL=0.8
% IC: Square Wave
% BC: Periodic
% CFL=0.8 
% dx=constant
% algorithm = Lax-Wendroff
clear all; close all; clc
SF=1;
range=[0,pi];
imax=401;
xpeak=0;
xc=NaN;
BC='Periodic';
IC='Square';
tfinal=imax;
nmax=1000;
cfl=0.8;

[index,x,dx,u,t]=Runner(range,imax,xpeak,xc,SF,BC,IC,tfinal,nmax,cfl);
%Create figure
figure1 = figure('Name','CP3.3_FigurePartE_SqWaveCFL0.8','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
xlim(axes1,[0 max(x)]);
ylim(axes1,[-.3 1.3]);
box(axes1,'on');
hold(axes1,'on');
xlabel(axes1,'x');
ylabel(axes1,'u(x)');
legend(axes1,'u(0,x)','u(imax\Deltat/4,x)')
title(axes1,{'Square Wave Propogation of u(t,x) with Periodic BC and Lax-Wendroff Algorithm (CFL=0.8)'});
plot1 = plot(x,u(1,:),x,u(101,:),'Parent',axes1);
set(plot1(1),'DisplayName','u(0,x)');
set(plot1(2),'DisplayName','u(imax\Deltat/4,0)');
legend(axes1,'show');
%% (E) Square Wave @ CFL=0.6
%IC: Square Wave
%BC: Periodic
%CFL=0.6
%dx=constant
% xc=N/A; BC=du/dx=0 at x=min,xmax
% algorithm = Lax-Wendroff
clear all; close all; clc
SF=1;
range=[0,pi];
imax=401;
xpeak=0; 
xc=NaN; % only needed for  SF=/=1
BC='Periodic';
IC='Square';
tfinal=imax;
nmax=1000;
cfl=0.6;

[index,x,dx,u,t]=Runner(range,imax,xpeak,xc,SF,BC,IC,tfinal,nmax,cfl);
%Create figure
figure1 = figure('Name','CP3.3_FigurePartE_SqWaveCFL0.6','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
xlim(axes1,[0 max(x)]);
ylim(axes1,[-.3 1.3]);
box(axes1,'on');
hold(axes1,'on');
xlabel(axes1,'x');
ylabel(axes1,'u(x)');
title(axes1,{'Square Wave Propogation of u(t,x) with Periodic BC and Lax-Wendroff Algorithm (CFL=0.6)'});
plot1 = plot(x,u(1,:),x,u(101,:),'Parent',axes1);
set(plot1(1),'DisplayName','u(0,x)');
set(plot1(2),'DisplayName','u(imax\Deltat/4,0)');
legend(axes1,'show');
%% (E) Gauss Wave @ CFL=0.8
%IC: Gauss Wave
%BC: Periodic
%CFL=0.8 
%dx=constant
%algorithm = Lax-Wendroff
clear all; close all; clc
SF=1;
range=[0,pi];
imax=401;
xpeak=0.5; %center of Gauss Functino
xc=NaN;
BC='Periodic';
IC='Gauss';
tfinal=imax;
nmax=1000;
cfl=0.8;

[index,x,dx,u,t]=Runner(range,imax,xpeak,xc,SF,BC,IC,tfinal,nmax,cfl);
%Create figure
figure1 = figure('Name','CP3.3_FigurePartE_GaussWaveCFL0.8','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
xlim(axes1,[0 max(x)]);
ylim(axes1,[0.7 2]);
box(axes1,'on');
hold(axes1,'on');
xlabel(axes1,'x');
ylabel(axes1,'u(x)');
title(axes1,{'Gauss Wave Propogation of u(t,x) with Periodic BC and Lax-Wendroff Algorithm (CFL=0.8)'});
plot1 = plot(x,u(1,:),x,u(101,:),'Parent',axes1);
set(plot1(1),'DisplayName','u(0,x)');
set(plot1(2),'DisplayName','u(imax\Deltat/4,0)');
legend(axes1,'show');
%% (E) Gauss Wave @ CFL=0.6
%IC: Gauss Wave
%BC: Periodic
%CFL=0.6
%dx=constant
%Gaussian function centered at x = 0.5 using a uniform grid 
clear all; close all; clc
SF=1;
range=[0,pi];
imax=401;
xpeak=0.5;
xc=NaN;
BC='Periodic';
IC='Gauss';
tfinal=imax;
nmax=1000;
cfl=0.6;

[index,x,dx,u,t]=Runner(range,imax,xpeak,xc,SF,BC,IC,tfinal,nmax,cfl);
%Create figure
figure1 = figure('Name','CP3.3_FigurePartE_GaussWaveCFL0.6','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
xlim(axes1,[0 max(x)]);
ylim(axes1,[0.7 2]);
box(axes1,'on');
hold(axes1,'on');
xlabel(axes1,'x');
ylabel(axes1,'u(x)');
title(axes1,{'Gauss Wave Propogation of u(t,x) with Periodic BC and Lax-Wendroff Algorithm (CFL=0.6)'});
plot1 = plot(x,u(1,:),x,u(101,:),'Parent',axes1);
set(plot1(1),'DisplayName','u(0,x)');
set(plot1(2),'DisplayName','u(imax\Deltat/4,0)');
legend(axes1,'show');
