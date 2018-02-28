% 30.10.17. Demo of simulating a time-continuous dynamical system and
% taking a Poincare section of it in phase space

clear 
close all

% Parameters of L84
a = 0.25;
b = 4;
F = 8; % summer - 6, winter - 8
G = 1;
Gp = 0;

% Random initial condition
x0 = [rand 0 0];

% Function handle
odefun = @(t,x)L84_rhs(t,x,a,b,F,G,Gp);

% Poincare section
eventfun = @(t,x)Poincare(t,x); 
options = odeset('Events',eventfun,'RelTol',1e-6);
tic
[T,X,TE,XE,IE] = ode45(odefun,[0 1e3],x0,options);
                               % 1e3 takes 3 sec
toc

% Simple plot
figure
plot(XE(:,1),XE(:,2),'.','MarkerSize',6)
xlabel('x'); ylabel('y');