function [m,h] = Optimal_Perturbation(M)
n=length(M);

%Step 1
[V,D] = eigs(M,1);
h = V;
h = h/sum(h);

%Step 2
B = triu(ones(n))-diag([1:n-1],-1);
B(:,n) = [];
B = sparse(normc(B));

%Step 3
n1 = length(find(M==0));
U = zeros(n,n^2-(n+n1));
j1 = 1;
j2 = 0;
for i=1:n
    R = find(M(:,i)==0);
    r = length(R);
    if r~= n-1
        B_i = zeros(n,n-r-1);
        R2 = setdiff([1:n],R);
        r2 = length(R2);
        B_i(R2,:) = B(1:r2,1:(r2-1));
        j2 = j2+n-r-1;
        U(:,j1:j2) = h(i)*B_i;
        j1 = j2+1;
    end
end

M_inf = h*ones(1,n);
Q = inv(eye(n)-M+M_inf);
U = Q*U;

%Step 4
[U2,D2,V2] = svds(U,1);

%Step 5
m = sparse(n,n);
j1=1;
j2=0;
for i=1:n
    R = find(M(:,i)==0);
    r = length(R);
    j2=n-r-1+j2;
    if r~= n-1
        B_i = zeros(n,n-r-1);
        R2 = setdiff([1:n],R);
        r2 = length(R2);
        B_i(R2,:) = B(1:r2,1:(r2-1));
        m(:,i) = B_i*V2(j1:j2);
    else
        m(:,i) = sparse(n,1);
    end
end