% 31.10.17. Computation of the maximal Lyapunov exponent from the
% separation of nearby trajectories.

clear 
close all

% Parameters of L84
a = 0.25;
b = 4;
F = 8; % summer - 6, winter - 8
G = 1;
Gp = 0;

% Switch between integrators; 1 for Matlab's ode45, 2 for user defined
% oderk4 and rk4 
sw = 1;

x0 = [rand 0 0];

% Function handle
odefun = @(t,x)L84_rhs(t,x,a,b,F,G,Gp);
options = odeset('RelTol',1e-10);

% Simulation of a long trajectory
dt = 1e-2; % time increment for output
n = 1e5; % number of time steps for output
if sw
    [T,X] = ode45(odefun,[0:dt:dt*n],x0,options);
else
    [T,X] = oderk4(odefun,dt,n,x0);
end

% Simulation of many short trajectories; ICs are chosen randomly from the
% long trajectory
nic = 100; % number of initial conditions
x0 = X(randsample(n,nic),:)'; % mind the transposition
x0 = x0(:);
n = 2e4; % number of time steps for output
if sw
    [T,X1] = ode45(odefun,[0:dt:dt*n],x0,options);
else
    [T,X1] = oderk4(odefun,dt,n,x0);
end

% Perturb each IC by a tiny number, close to machine precision
x0 = x0 + 1e-10*rand(size(x0));
if sw
    [T,X2] = ode45(odefun,[0:dt:dt*n],x0,options);
else
    [T,X2] = oderk4(odefun,dt,n,x0);
end

% Obtain the distance of the reference and pertubed trajectories
d = sqrt((X1(:,      1:  nic) - X2(:,      1:  nic)).^2 + ...
         (X1(:,  nic+1:2*nic) - X2(:,  nic+1:2*nic)).^2 + ...
         (X1(:,2*nic+1:end  ) - X2(:,2*nic+1:end  )).^2);

ord  = mean(log(d),2);
absc = T;



figure;  plot(log(d))

hold on; plot(ord,'LineWidth',4,'Color','k')

% Fitting straight line
rg = 3e3:1e4; % fit range 
[p,S] = polyfit(absc(rg),ord(rg),1);
%title(['MLE = ' num2str(p(1))])
%ylabel('log(distance)'); xlabel('time');