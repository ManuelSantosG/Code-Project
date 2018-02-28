function y = mixexpcdf(x,p,mu1,mu2)


y=(1-p)*expcdf(x,mu1)+p*expcdf(x,mu2);


end

