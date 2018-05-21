s = 10; rh = 28; b = 8/3;                % the Lorenz system
epsilon=0.1;
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3)+epsilon*x(:,1)...
          x(:,1).*x(:,2)-b*x(:,3)];
h = 0.001; n = 1; f = @(x) rk4(v,x,h,n); % f is the flow


x_i=[1 1 1];
fold=x_i;
NN=10000000;
obs1=x_i(1)^2;
obs2=x_i(2)^2;
obs3=x_i(3)^2;
obs4=x_i(3);
for i = 1 :NN-1;
    fnew=f(fold);
    obs1=obs1 + fnew(1)^2;
    obs2=obs2 + fnew(2)^2;
    obs3=obs3 + fnew(3)^2;
    obs4=obs4 + fnew(3);
    fold=fnew;
end
obs1=obs1/NN
obs2=obs2/NN
obs3=obs3/NN
obs4=obs4/NN