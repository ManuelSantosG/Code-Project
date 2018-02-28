% 30.10.2017. Demonstration for the computation of Lyapunov exponents of
% time-continuous dynamical system.

close all
clear

% Parameters (L84)
a = 0.25;
b = 4;
F = 8;
G = 1;
Gp = 0;

% Function handle for the ODE
rhs_ode_var = @(t,x)L84_var(t,x,a,b,F,G,Gp);

% Initial conditions
x0 = [rand 0 0];

t0 = 0; % initial time
h_rn = 10; 
    % It must be small enough. h_rn=1 seems to balance well accuracy and
    % numerical efficiency.
tf = 5000; % final time

% Number crunching
tic
[T,Res] = lyapunov(length(x0),rhs_ode_var,@ode45,t0,h_rn,tf,x0,0);
toc

% Visualise results: finite-time approximations of the LE's
figure;
plot(T,Res);
xlabel('Time'); ylabel('Lyapunov exponents');

% Calculation of the Kaplan-Yorke dimension DKY
lambdas = sort(Res(end,:),'descend'); 
    % This includes 0 LE(s) too.
csls = cumsum(lambdas);
k = find(csls<0,1)-1; 

DKY = k + csls(k)/abs(lambdas(k+1));
    % lambdas(k+1) will never be 0
title(['Kaplan-Yorke dimension = ' num2str(DKY)]);