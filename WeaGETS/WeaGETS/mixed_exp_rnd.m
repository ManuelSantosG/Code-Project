function [X]=mixed_exp_rnd(value,alpha,beta1,beta2)

% [X]=mixed_exp_rnd(value,n,alpha,beta1,beta2);
%
% This function computes a random sample of size 'n' from the three parameters
% describing a mixed-exponential distributions
%
% Input
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'n' : integer number = length of the time series
% 'alpha' : real number (for choosing one of both means) 
% 'beta1' : real number (= mean of first exponential distribution)
% 'beta2' : real number (= mean of second exponential distribution)
%
% Output
% 'X' : vector of length 'n'
% 
% Vincent MORON
% december 2004

P=rand(1);
if alpha>value  
    X=expinv(P,beta1);  %exponential dist
else 
    X=expinv(P,beta2);  %exponential dist
end


