function MatrixPowers(A);

n=50;
l=zeros(1,n);
[n1,n2]=size(A);
M=eye(n1);
for i=1:n;
    M=M*A;
    l(i)=norm(M);
end

plot(l)
end