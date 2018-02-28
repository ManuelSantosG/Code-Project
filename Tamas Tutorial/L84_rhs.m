function xprime = L84_rhs(t,x,a,b,F,G,Gp)

% Vectorised right hand side of ODE, suitable for handling ensemble 

n = length(x)/3;

xprime = [-x(n+1:2*n).^2 - x(2*n+1:end).^2 - a*x(1:n) + a*F;
    x(1:n).*x(n+1:2*n) - b*x(1:n).*x(2*n+1:end) - x(n+1:2*n) + G;
    x(1:n).*x(2*n+1:end) + b*x(1:n).*x(n+1:2*n) - x(2*n+1:end) + Gp];