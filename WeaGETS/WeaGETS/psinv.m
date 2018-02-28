function [rain_values]=psinv(R,par1,par2,par3)
% this program is to calculate the rainfall value with the Pearson III
% distribution 

rain_values=par1+(2*par2/par3)*(((par3/6)*(R-(par3/6))+1)^3-1);

end
