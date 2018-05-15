function [l]=MatrixPowers(A);
    tic

    n=100;
    l=zeros(2,n);
    [n1,n2]=size(A);
    M=eye(n1);
    for i=1:n;
        M=A^i;
        r=norm(M);
        l(1,i)=r;
        %l(1,i)=nthroot(r,i);
    end

    plot(l(1,:))
    hold on
    plot(l(2,:))
    n1;
    toc
end