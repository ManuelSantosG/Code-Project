%Ulam's method for the generator

F=@(x) sin(4*pi*x) + 1.1;
n=128;
Fx=F(0:1/n:1-1/n);
A=n*sparse(1:n,1:n,-Fx) + n*sparse(1:n,[2:n 1],Fx);
[v,d]=eig(full(A'));
[d,I]=sort(diag(d));
f_n = @(x) abs(v(ceil(max(x*n,1)),I(1)));
L1=quadl(f_n,0,1);
x=0:1/n:1-1/n;
fff=[];
for i = 1:n
    fff=[fff f_n(x(i))];
end



 fplot(@(x)f_n(x)/L1,[0,1]);