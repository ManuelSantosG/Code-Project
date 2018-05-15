function [l]=MatrixExponential(A);
    tic
    n=500;
    T=50;
    t=linspace(0,T,n);
    [n1,n2]=size(A);
    M=eye(n1);
    l=zeros(1,n);
    for i=1:n;
        M=expm(A*t(i));
        r=norm(M,1);
        l(i)=r;
        %l(1,i)=nthroot(r,i);
    end

    plot(l);
    l(end)
    
    toc
end