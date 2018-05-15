
function psi1=Response_Operator(M,m);
    n=100;
    l=size(M,2);
    psi1=M;
    Mn=eye(l);
    for i = 1 : n;
        Mn=Mn*M;
        psi1=psi1 + Mn;
    end
    psi1=psi1';
    psi1=m'*psi1;
    
end