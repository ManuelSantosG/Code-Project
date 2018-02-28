function [ParHat LogLik] = mixexpfit(x)
% parameter estimates for mixed exponential distribution by MLE
%##########################################################################
%
% mixexpfit estimate the parameters of mixed exponential distribution by
% maximum likelihood estimation method.
%
% INPUTS:
% x -> iid random sample (row or column vector or matrix)
%
% OUTPUTS:
% ParHat -> estimated parameters vector [mixing factor, mu1, mu2]
% LogLik -> loglikelihood function evaluated at ParHat
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
% See also: mixexprnd
%


% input check
%--------------------------------------------------------------------------
error(nargchk(1,1,nargin));

x=x(:);
if any(isnan(x))
    error('NaN existed in x.');
end


% assign initial parameters randomly
%--------------------------------------------------------------------------
InitPar=[log(rand/(1-rand)) 10*rand 20*rand];


% define an inline function
%--------------------------------------------------------------------------
objfun=@(Par,x) -sum(log(mixexppdf(Par,x)));


% search parameters minimize the negative log-likelihood function
%--------------------------------------------------------------------------
options=optimset('MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);
p=fminsearch(objfun,InitPar,options,x);
ParHat=[InvMixFactor(p(1)) p(2:3)];
LogLik=-objfun(p,x);

% Note: since mixing factor must be a real value falling between 0 and 1,
% to avoid the bounded optimization problem, the mixing factor was
% transformed by log(MixFactor/1-MixFactor). That is why InvMixFactor is
% necessary.


%% subroutine 1
    function y=mixexppdf(Par,x)
        % compute the probability density of mixed exponential distribution
        
        p=InvMixFactor(Par(1));
        mu1=Par(2);
        mu2=Par(3);
        
        y=p*exppdf(x,mu1)+(1-p)*exppdf(x,mu2);
        
    end


%% subroutine 2
    function mixp=InvMixFactor(TransP)
        % backtransform the tranformed mixing factor into original scale
        % TransP=log(mixp/1-mixp)
        
        mixp=exp(TransP)/(1+exp(TransP));
        
    end


end

