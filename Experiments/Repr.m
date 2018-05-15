
function psi1=Response_Operator(M,m);
    n=50;
    l=size(M,2);
    psi1=eye(l);
    for i = 1 : 50;
        psi1=psi1*M;
    end
    psi1=psi1';
    psi1=m'*psi1;