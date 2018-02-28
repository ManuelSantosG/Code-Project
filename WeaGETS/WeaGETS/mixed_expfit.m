function [best_estimates]=mixed_expfit(x,maxiter);

% [best_estimates]=mixed_expfit(x,maxiter);
%
% Computation of the parameters [alpha, beta1 and beta2] of a mixed 
% exponential distribution with 3 parameters (one for the the mixing
% probability and two for the exponential ditribution). 
%
% Input
% 'x' : vector of real number describing the sample of data
% 'maxiter' : scalar integer giving the maximum number of iterations
%
% Output
% 'best_estimates' : vector of three real number giving the estimates of a
% mixed exponential PDF fitting the values of 'x'.
%
% program coded by Tom Lane (tom.lane@mathworks.com)
% March 2005

rand('seed',sum(100*clock));
randn('seed',sum(100*clock));

%disp(sprintf('%d iterations with random starting values',maxiter))
fmax = -inf;
for j=1:20
    pstart = rand;
    start = [log(pstart/(1-pstart)), 10*rand, 20*rand];
    p = fminsearch(@objfun,start,[],x);
    fval = -objfun(p,x);
%    if (j<3 | j>maxiter-2)
%       disp([j,fval,logit(p(1)),p(2),p(3)])
%    elseif j==3
%       disp('    ...')
%    end
    if fval>fmax
        fmax = fval;
        pbest = p;
    end
end

%best_log_likelihood = fmax
best_estimates = [logit(pbest(1)),pbest(2),pbest(3)];

function f = objfun(param,x)
% objective function is the negative log likelihood
f = -sum(log(mixdens(param,x)));

function d = mixdens(param,x)
% compute the mixture density
% param(1) is log(p/(1-p)) where p is the proportion in the original
% mixture formulation, so param(1) does not need to be bounded
p = logit(param(1));
a = param(2);  % mean of 1st exponential component
b = param(3);  % mean of 2nd exponential component
d = p*exppdf(x,a) + (1-p)*exppdf(x,b);

function l=logit(p)
% undo the log(p/(1-p)) transformation
l = exp(p)/(1+exp(p));