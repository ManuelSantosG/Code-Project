% 31.10.17. Demonstration of computing the box dimension of a
% fractal set living in the plane, e.g., the trajectory of a discrete time
% dynamical system, i.e., a map.

clear
close all

load XE

xsl = XE(:,1);
ysl = XE(:,2);

figure
plot(XE(:,1),XE(:,2),'.','MarkerSize',2)
xlabel('x'); ylabel('y');

eps = 2.^-linspace(1,10,100); % resolution/size of box
N = 0*eps; % Number of boxes to fully cover attractor

for i = 1:length(eps)
    % Which 'x and y slots' is each point in?
    xb = floor(xsl/eps(i));
    yb = floor(ysl/eps(i));

    % Now we have to eliminate redundant pairs of (xb,yb)
    nx = max(xb) - min(xb) + 1; % number of x slots that span the object

    % Identify a box with a single number using the concept of digits
    xyb = xb + nx*yb;
    
    N(i) = length(find(diff(sort(xyb))));
end

ord = log2(N);
absc = log2(1./eps);

figure
plot(absc,ord,'o','Color','magenta')
xlabel('log_2(1/\epsilon)')
ylabel('log_2(N(\epsilon))')

% Fitting straight line
rg = 20:60; % fit range 

[p,S] = polyfit(absc(rg),ord(rg),1);
D0 = p(1);
f = polyval(p,absc(rg));
RMSE = sqrt(sum((f - ord(rg)).^2)/length(f));
title(['D_0 = ' num2str(D0) ', RMSE = ' num2str(RMSE)])
hold on;
plot(absc,polyval(p,absc),'Color','k')
plot(absc(rg),f,'Color','g','LineWidth',2) % range where the fit took place