function [value,isterminal,direction] = Poincare(t,x)

n = length(x)/3;

value = x(2*n+1:end);
isterminal = 0*value;
direction = ones(size(value));