function [t,x] = oderk4(odefun,h,n,x0)

t = (0:h:h*n)';
nt = length(t);
x = zeros(length(x0),nt);
x(:,1) = x0;

hand = waitbar(0,'Simulation in progress...');
for i = 2:nt
    x(:,i) = RK4(odefun,t(i),h,x(:,i-1));
    waitbar(i/nt)
end
x = x';
close(hand) 
