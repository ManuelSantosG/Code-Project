s = 10; rh = 28; b = 8/3;                % the Lorenz system
epsilon=0.1;
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3)+epsilon*x(:,1)...
          x(:,1).*x(:,2)-b*x(:,3)];
h = 0.01; n = 10; f = @(x) rk4(v,x,h,n); % f is the flow


x_i=[1 1 1];
fold=x_i;
NN=1000;
CC=x_i(1)^2;
for i = 1 :NN-1;
    fnew=f(fold);
    CC=CC + fnew(1)^2;
    fold=fnew;
end
CC=CC/NN