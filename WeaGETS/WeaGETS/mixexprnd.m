function R = mixexprnd(p,mu1,mu2,n)
% mixed exponential distribution random numbers
%##########################################################################
%
% generates random numbers from mixed exponential distribution
%
% INPUTS:
% p -> mixing factor
% mu1 -> parameter for the 1st exponential distribution
% mu2 -> parameter for the 2nd exponential distrbution
% n -> the number of random variates to be generated
%
% OUTPUT:
% R -> column vector of the generated random variates
%
% Reference:
% D.A. Woolhiser, J. Roldan, Stochastic Daily Precipitation Models 2: a
% Comparison of Distribution of Amounts. Water Resources Research, 1982.
%
%
%
% This function serves for StoRainGen toolbox.
% 
% Latest Updated: 01/05/2011
%
% Copyright © 01/01/2011 Lichao, BAEN, TAMU
%
%##########################################################################
%
% See also: mixexpfit
%

% input check
%--------------------------------------------------------------------------
error(nargchk(3,4,nargin));

if nargin < 4
    n=1;
end

% define an inline function
%--------------------------------------------------------------------------
objfun=@(x,Par) Par(1)-1+Par(2)*exp(-x/Par(3))+(1-Par(2))*exp(-x/Par(4));


% find the zero of the mixture exponential CDF
%--------------------------------------------------------------------------
options=optimset('MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
R=zeros(n,1);
U=rand(n,1);
for i=1:n
    Par=[U(i) p mu1 mu2];
    R(i)=fzero(objfun,rand,options,Par); % using random initial quantile
end


end

