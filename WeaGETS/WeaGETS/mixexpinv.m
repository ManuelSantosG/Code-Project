function x = mixexpinv(Px,p,mu1,mu2)

% define an inline function
%--------------------------------------------------------------------------
objfun=@(x,Par) Par(1)-1+Par(2)*exp(-x/Par(3))+(1-Par(2))*exp(-x/Par(4));

% find the zero of the mixture exponential CDF
%--------------------------------------------------------------------------
Px=Px(:);
n=length(Px);
options=optimset('MaxIter',10000,'TolFun',1e-10,'TolX',1e-10);
x=zeros(n,1);
for i=1:n
    Par=[Px(i) p mu1 mu2];
    x(i)=fzero(objfun,rand,options,Par); % using random initial quantile
end

end

