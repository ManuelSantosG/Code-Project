% 30.10.17. Demo of simulating a time-continuous dynamical system

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

% Simulation
dt = 1e-2; % time increment for output
n = 2^15; % number of time steps for output; some power of 2 
options = odeset('RelTol',1e-10);
[T,X] = ode45(odefun,[0:dt:dt*n],x0,options);

% Visualisation of results
% Trajectory in phase space
figure
plot3(X(:,1),X(:,2),X(:,3))
xlabel('x'); ylabel('y'); zlabel('z'); 
% Time series of components
figure
subplot(3,1,1); plot(T,X(:,1)); ylabel('x'); 
subplot(3,1,2); plot(T,X(:,2)); ylabel('y'); 
subplot(3,1,3); plot(T,X(:,3)); ylabel('z'); xlabel('time');

% Power spectrum 
% https://uk.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html 
fs = 1/dt; 
xdft2 = fft(X(:,1));
xdft2 = xdft2(1:n/2+1);
psdx2 = (1/(fs*n)) * abs(xdft2).^2;
psdx2(2:end-1) = 2*psdx2(2:end-1);
f2 = 0:fs/n:fs/2;

[psdx3, f3] = periodogram(X(:,1),[], n, fs);
figure; semilogx(f2,10*log10(psdx2),f3,10*log10(psdx3)); 
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')