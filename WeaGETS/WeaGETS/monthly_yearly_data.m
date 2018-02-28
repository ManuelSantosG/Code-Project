function monthly_yearly_data(filenamein)
% this program is to calculate the monthly and yearly precipitation and 
% yearly averaged Tmax and Tmin for the observed data

load(filenamein) % load the observed data

% set missing values (-999) with NaN

nn=find(P==-999);
P(nn)=NaN;
nn=find(Tmax==-999);
Tmax(nn)=NaN;
nn=find(Tmin==-999);
Tmin(nn)=NaN;

n=size(P,1);

% calculate the observed monthly precip
monthly_observed_P=zeros(n,12);
for i=1:n
    monthly_observed_P(i,1)=nansum(P(i,1:31));
    monthly_observed_P(i,2)=nansum(P(i,32:59));
    monthly_observed_P(i,3)=nansum(P(i,60:90));
    monthly_observed_P(i,4)=nansum(P(i,91:120));
    monthly_observed_P(i,5)=nansum(P(i,121:151));
    monthly_observed_P(i,6)=nansum(P(i,152:181));
    monthly_observed_P(i,7)=nansum(P(i,182:212));
    monthly_observed_P(i,8)=nansum(P(i,213:243));
    monthly_observed_P(i,9)=nansum(P(i,244:273));
    monthly_observed_P(i,10)=nansum(P(i,274:304));
    monthly_observed_P(i,11)=nansum(P(i,305:334));
    monthly_observed_P(i,12)=nansum(P(i,335:365));
end
save('monthly_observed_P','monthly_observed_P')

% calculate the observed yearly precip,Tmax and Tmin
yearly_observed_P=zeros(n,1);
for i=1:n
    yearly_observed_P(i,1)=nansum(P(i,:));
    yearly_observed_tmax(i,1)=nanmean(Tmax(i,:));
    yearly_observed_tmin(i,1)=nanmean(Tmin(i,:));
end
yearly_observed_P=yearly_observed_P';
yearly_observed_tmax=yearly_observed_tmax';
yearly_observed_tmin=yearly_observed_tmin';
save('yearly_observed_P','yearly_observed_P')
save('yearly_observed_tmax','yearly_observed_tmax')
save('yearly_observed_tmin','yearly_observed_tmin')

