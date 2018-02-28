function y = mixexppdf(x,p,mu1,mu2)


y=(1-p)*exppdf(x,mu1)+p*exppdf(x,mu2);

end

