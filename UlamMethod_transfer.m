F = @(x) sin(4*pi*x)+1.1;
n = 128; x = 1/(2*n^2):1/(n^2):1-1/(2*n^2);
[t,y] = ode23(@(t,x) F(x),[0 2/n],x');
I = ceil(max(n*mod(y(end,:),1),1));
J = reshape(ones(n,1)*(1:n),1,n*n);
P = sparse(I,J,1/n,n,n);
[v,d] = eig(full(P));
f_n = @(x) abs(v(ceil(max(x*n,1)),1));
L1 = quadl(f_n,0,1); fplot(@(x)f_n(x)/L1,[0,1])
axis([-1 1 -1 1])