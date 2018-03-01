function [idistr,parp000,parp010,parp100,parp110,parlbd,parnu,A,B,aC0,aC1,aC2,...
    aD1,aD2,sC0,sC1,sC2,sD1,sD2,PrecipThreshold,MarkovChainOrder]=analyzer_order2_norad...
    (filenamein,PrecipThreshold,idistr,harm_order)
%
% THIS PROGRAM ANALYSES DAILY TIME SERIES OF PRECIPITATION AND MAX AND MIN
% TEMPERATURES 

% the preciptiation occurrence is generated using the second order markov
% chain and the parameters are smoothed from 14 days to daily scale with
% this weather generator

% monthly matrices of daily precip are read in, starting in January
% feb 29 omitted. 
%
% The following matrices must be included in the matlab file named
% 'filenamein' (string):
%
% P = daily precipitation in mm in the matrix format [nyears x 365]
% Tmax = maximal daily temperature in C in the matrix format [nyears x 365]
% Tmin = minimal daily temperature in C in the matrix format [nyears x 365]
% yearT =  vector of years with Tdata
% yearP = vector of years with Precip data
% missing values must be -999
% 
% 'filenameout' is a string for the filename to which results will be
% written to (.mat).  This filename can be used by the generator for
% generation of synthetic data
% 'station_name' is a string for the station name for graphing purposes
%
% PrecipThreshold: amount of precipitation that is considered to determine 
% whether a given day is wet
% idistr: 1=exponential, 2=gamma  (distribution function to be fitted to precip) 
% harm_order: order of harmonics used to fit the yearly distribution of the
% markov parameters and distribution functions.

%Indicate wich order was used to analyse the data.
%This is necessary to save this attributes in the output file.
%It will determine the generator order to be used.
MarkovChainOrder = 2;

h = waitbar(0,'Data analysis...');

load(filenamein);
%
% we start the analysis of the precipitation time series
%
% generate matrix of precip - no precip (PnP)

[Y,D]=size(P);      % years, days
PnP=zeros(Y,D);
[k]=find(P>PrecipThreshold);        % finds days with significant precipitation
PnP(k)=1;           % a value of 1 means that significant precipitation occured
%
% adjust precipitation to threshold: substract threshold to precipitation
% values so to be able to use the 1-parameter exponential funtion of
% 2-parameter gamma funtion.  The threshold will be added back by the
% generator

[k2]=find(P<=PrecipThreshold & P>0); 
P(k2)=0;
P(k)=P(k)-PrecipThreshold;
%
% put NaN for missing values instead of -999
%
[kn]=find(P<0);
PnP(kn)=NaN;
%
%  calculate transition state a00, a01, a10, a11 (Woolhiser and Pegram, 1979)
%
% matrix 'a' is a [4 365] matrix that contains the number of years over Y
% where a00 (line 1) a01 (line 2) a10 (line 3) and a11 (line 4) was
% observed.  As such, the sum of each column of matrix 'a' is equal to Y
% (in fact it is equal to Y-Days_without_data)

[a] = transition2(PnP,Y,D);

waitbar(0.05,h);
%
% calculate average p00 and p10 using maximum likelihood estimator
% on 14-days periods (Woolhiser and Pegram, 1979)
%
n14=round(D/14);
for i=1:n14
   A000(i)=sum(a(1,14*(i-1)+1:14*i));
   A001(i)=sum(a(2,14*(i-1)+1:14*i));
 
   A010(i)=sum(a(3,14*(i-1)+1:14*i));
   A011(i)=sum(a(4,14*(i-1)+1:14*i));
  
   A100(i)=sum(a(5,14*(i-1)+1:14*i));
   A101(i)=sum(a(6,14*(i-1)+1:14*i));
 
   A110(i)=sum(a(7,14*(i-1)+1:14*i));
   A111(i)=sum(a(8,14*(i-1)+1:14*i));
end

ap000=A000./(A001+A000);
ap001=1-ap000;

ap010=A010./(A010+A011);
ap011=1-ap010;

ap100=A100./(A101+A100);
ap101=1-ap100;

ap110=A110./(A110+A111);
ap111=1-ap110;

waitbar(0.1,h);

%
% Fourier analysis of daily values of p00 an p10 using a maximum likelihood approach
% according to Woolhiser and Pegram (1979)
%
% we first estimate the parameters of the Fourier series using a least square fit with
% the ap00 and ap10 values as a starting point for the maximum likelihood approach
%
% a maximum of four harmonics is considered accurate to model the variations of p00 and p10
%
t=1:n14;
t=t';
T14=n14/(2*pi);
%
% first, p000
%
[initphi(1,:)] = fourier14(ap000,t,T14);
%
% first, p010
%
[initphi(2,:)] = fourier14(ap010,t,T14);
%
% then, p100
%
[initphi(3,:)] = fourier14(ap100,t,T14);
%
% first, p110
%
[initphi(4,:)] = fourier14(ap110,t,T14);

waitbar(0.2,h);
%
% an optimisation is performed to obtain the Fourier parameters using
% a maximum likelihood approach with the daily time series
%
% no harmonics
%
suma000=sum(a(1,:));
suma001=sum(a(2,:));

suma010=sum(a(3,:));
suma011=sum(a(4,:));

suma100=sum(a(5,:));
suma101=sum(a(6,:));

suma110=sum(a(7,:));
suma111=sum(a(8,:));

C000opt=suma000/(suma000+suma001);

C010opt=suma010/(suma010+suma011);

C100opt=suma100/(suma100+suma101);

C110opt=suma110/(suma110+suma111);

optphi=initphi;
optphi(1,1)=C000opt;

optphi(2,1)=C010opt;

optphi(3,1)=C100opt;

optphi(4,1)=C110opt;
%
%  likelihood function with zero harmonic
%
U0(1)=log(C000opt)*suma000+log(1-C000opt)*suma001;

U1(1)=log(C010opt)*suma010+log(1-C010opt)*suma011;

U2(1)=log(C100opt)*suma100+log(1-C100opt)*suma101;

U3(1)=log(C110opt)*suma110+log(1-C110opt)*suma111;

% optimizing all 4 transition probabilities.  Proceed one harmonic at a time
%
n=1:D;
T=D/(2*pi);
%
% optimizing p000
%
for i=1:4
    pt=initphi(1,2*i:2*i+1);
    if i==1
        funk=@mle_p1;
    elseif i==2
        funk=@mle_p2;
    elseif i==3
        funk=@mle_p3;
    else
        funk=@mle_p4;
    end
    a1=a(1,:);
    a2=a(2,:);
    [ptopt,U0(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        optphi(1,1:2*i-1),a1,a2,n,T);
    optphi(1,2*i:2*i+1)=ptopt;
end
waitbar(0.275,h);

%
% optimizing p010
%
for i=1:4
    pt=initphi(2,2*i:2*i+1);
    if i==1
        funk=@mle_p1;
    elseif i==2
        funk=@mle_p2;
    elseif i==3
        funk=@mle_p3;
    else
        funk=@mle_p4;
    end
    a3=a(3,:);
    a4=a(4,:);
    [ptopt,U1(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        optphi(2,1:2*i-1),a3,a4,n,T);
    optphi(2,2*i:2*i+1)=ptopt;
end
waitbar(0.425,h);
%
% optimizing p100
%
for i=1:4
    pt=initphi(3,2*i:2*i+1);
    if i==1
        funk=@mle_p1;
    elseif i==2
        funk=@mle_p2;
    elseif i==3
        funk=@mle_p3;
    else
        funk=@mle_p4;
    end
    a5=a(5,:);
    a6=a(6,:);
    [ptopt,U2(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        optphi(3,1:2*i-1),a5,a6,n,T);
    optphi(3,2*i:2*i+1)=ptopt;
end
waitbar(0.575,h);
%
% optimizing p110
%

for i=1:4
    pt=initphi(4,2*i:2*i+1);
    if i==1
        funk=@mle_p1;
    elseif i==2
        funk=@mle_p2;
    elseif i==3
        funk=@mle_p3;
    else
        funk=@mle_p4;
    end
    a7=a(7,:);
    a8=a(8,:);
    [ptopt,U3(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        optphi(4,1:2*i-1),a7,a8,n,T);
    optphi(4,2*i:2*i+1)=ptopt;
end
waitbar(0.6,h);

%  we now process the precipitation amounts
%  an exponential distribution (Richardson, 1981), or a 2-parameter Gamma
%  distribution is used
%
%  we first create a precipitation matrice containing positive
% precipitation amounts when precipitation has been recorded and
% zero values (when there is no precipitation or when data is missing)
% otherwise.  This requires to replace NaN values in PnP matrice
% by 0 values.  Position of NaN values in vector kn
%
PnP(kn)=0;
PP=P.*PnP;
sPtot=sum(PP); %total precipitation amounts on given calendar day for all years analysed
stot=sum(PnP); % # of days with precip. on a given calendar day for all years analysed
sPtot2=sum(PP.^2);
%
%  select exponential or gamma distribution
%
if idistr == 1   % exponential distribution
    %
    %  mle of parameter lambda of the exponential distribution taken over 14-day periods
    %
    for i=1:n14
        alambda(i)=sum(stot(14*(i-1)+1:14*i))/sum(sPtot(14*(i-1)+1:14*i));
    end
    %
    % Fourier analysis on 14-day lambda values uding least square fit approach
    %  maximum of 4 harmonics
    %
    [phi2] = fourier14(alambda,t,T14);
    %
    % optimizing lambda using maximum likelihood approach
    %
    % first, 0 harmonics
    %
    optphi2=zeros(1,length(phi2));
    optphi2(1)=phi2(1);
    U8(1)=mle_lam0(optphi2,PP,PnP);
    U8(1)=-U8(1);
    %
    % then harmonics 1 to 4.  Nelder and Mead approach used
    %
    v1=ones(Y,1);
    for i=1:4
        pt=phi2(2*i:2*i+1);
        if i == 1
            funk=@mle_lam1;
        elseif i == 2
            funk=@mle_lam2;
        elseif i == 3
            funk=@mle_lam3;
        else
            funk=@mle_lam4;
        end
        [ptopt,U8(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        optphi2(1:2*i-1),v1,PP,PnP,n,T);
        optphi2(2*i:2*i+1)=ptopt;
    end
    parnu=-999;       % just so that the variable is defined to avoid errors in the SAVE portion

elseif idistr == 2   %  2-parameter gamma distribution
   
    kzero=find(PP==0);
    PP(kzero)=1e-20;
    logPP=log(PP);
    logPP(kzero)=0;
    %
    %  estimates of the lambda and nu parameters of the Gamma
    %  distribution taken over 14-day periods
    %
    slogPtot=sum(logPP);
    for i=1:n14
        aa(i)=sum(stot(14*(i-1)+1:14*i));
        xmean(i)=sum(sPtot(14*(i-1)+1:14*i))/aa(i);
        xvar(i)=(sum(sPtot2(14*(i-1)+1:14*i))-aa(i)*xmean(i)^2)/(aa(i)-1);
        meanlogx(i)=sum(slogPtot(14*(i-1)+1:14*i))/aa(i);
    end
    logxmean=log(xmean);
    y=logxmean-meanlogx;
    
    % MLM 
        nu=(1+sqrt(1+4*y/3))./(4*y);
        nu=(aa-3).*nu./aa;
   
    alambda=nu./xmean;
    [phi3]=fourier14(nu,t,T14);
    [phi4]=fourier14(alambda,t,T14);
    %
    % simultaneously optimizing lambda and nu by optimizing a mle function.
    %
    pt(1)=phi3(1);
    pt(2)=phi4(1);
    funk=@mle_gam0;
    [ptopt,U8(1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1e-5),...
        PP,logPP,PnP);
    pt(1)=ptopt(1);
    pt(2)=ptopt(2);
    v1=ones(Y,1);
    for i=1:4
        waitbar(0.6+0.075*i,h); 
        pt(4*i-1:4*i)=phi3(2*i:2*i+1);
        pt(4*i+1:4*i+2)=phi4(2*i:2*i+1);
        if i == 1
            funk=@mle_gam1;
        elseif i == 2
            funk=@mle_gam2;
        elseif i == 3
            funk=@mle_gam3;
        elseif i ==4
            funk=@mle_gam4;
        end
        [ptopt,U8(i+1),exitflag,output]=fminsearch(funk,pt,optimset('TolX',1) ...
            ,v1,PP,logPP,PnP,n,T);
      pt=ptopt;
        %ptopt
    end
    optphi3(1)=ptopt(1);
    optphi3(2:3)=ptopt(3:4);
    optphi3(4:5)=ptopt(7:8);
    optphi3(6:7)=ptopt(11:12);
    optphi3(8:9)=ptopt(15:16);
    optphi4(1)=ptopt(2);
    optphi4(2:3)=ptopt(5:6);
    optphi4(4:5)=ptopt(9:10);
    optphi4(6:7)=ptopt(13:14);
    optphi4(8:9)=ptopt(17:18);
end

   waitbar(0.96,h);
%
% results are plotted and incremental changes in U values are calculated
% to help select the number of significant harmonics
%
n=n';
if idistr==1
%
%  keep the number of significant
%  harmonics based on ham_order parameter    
%
    sh=[harm_order harm_order harm_order harm_order harm_order];
%
% store the resulting parameters in vectors for later use
%
    parp000=zeros(1,9);
    
    parp010=zeros(1,9);
   
    parp100=zeros(1,9);
   
    parp110=zeros(1,9);
    
    parlbd=zeros(1,9);
    
    % NB: keep fourier coefficients using least-square method for transition probabilities 
    
    parp000(1,1:1+2*sh(1))=initphi(1,1:1+2*sh(1));
    
    parp010(1,1:1+2*sh(2))=initphi(2,1:1+2*sh(2));
   
    parp100(1,1:1+2*sh(3))=initphi(3,1:1+2*sh(3));
   
    parp110(1,1:1+2*sh(4))=initphi(4,1:1+2*sh(4));
   
    parlbd(1:1+2*sh(5))=optphi2(1:1+2*sh(5));
else
    % gamma distribution
   
%  keep the number of significant
%  harmonics based on ham_order parameter
%    
    
    sh=[ harm_order harm_order harm_order harm_order harm_order harm_order];

%
% store the resulting parameters in vectors for later use
%
    parp000=zeros(1,9);
    
    parp010=zeros(1,9);
   
    parp100=zeros(1,9);
   
    parp110=zeros(1,9);
    parlbd=zeros(1,9);
    parnu=zeros(1,9);
    
    % NB: keep fourier coefficients using least-square method for
    % transition probabilities 
    
    parp000(1,1:1+2*sh(1))=initphi(1,1:1+2*sh(1));
    
    parp010(1,1:1+2*sh(2))=initphi(2,1:1+2*sh(2));
    
    parp100(1,1:1+2*sh(3))=initphi(3,1:1+2*sh(3));
   
    parp110(1,1:1+2*sh(4))=initphi(4,1:1+2*sh(4));
   
    parlbd(1:1+2*sh(5))=optphi4(1:1+2*sh(5));
    parnu(1:1+2*sh(6))=optphi3(1:1+2*sh(6));

end
%
% we start the analysis of the climatic (temperature, radiation) time series
%
%
%  PnP:  matrix with elements = 1 for wet days only (=0 for dry and missing
%  days)
%
%  PnPm1:  matrix with elements = 1 for dry days only (=0 when wet and
%  missing days).
%
PnPm1=1-PnP;  % here elements = 1 include dry and missing precipitation days
PnPm1(kn)=0;   % missing precipitation days are removed from the matrix
%
%  make sure that all matrices P and Tmin, Tmax and Rad are the same size
%
% take the years that are common only

[int_vec,ind_P,ind_T]=intersect(yearP,yearT);

P=P(ind_P,:);
PnP=PnP(ind_P,:);
PnPm1=PnPm1(ind_P,:);
Tmin=Tmin(ind_T,:);
Tmax=Tmax(ind_T,:);
[Y,D]=size(P);
%
% add 0.0001 to all temperatures and radiations equal to zero.
% this is required to efficiently and correctly
% calculate the mean and standard deviations
%
x=find(Tmax==0);
Tmax(x)=0.0001;
x=find(Tmin==0);
Tmin(x)=0.0001;
%
% separate the matrices Tmin, Tmax and Rad into wet and dry components
%
Tmaxw=Tmax.*PnP;
Tminw=Tmin.*PnP;
%
Tmaxd=Tmax.*PnPm1;
Tmind=Tmin.*PnPm1;
%
% calculate daily means and standard deviations for Tmaxw, Tmaxw, etc..
% missing values (-999) should not be included in the analysis.  A matrix
% of zeros and ones is created for each climate variable (Tmax, Tmin, Rad)
% to identify missing values.  Missing values are replaced by 0.
%
Mtmax=ones(Y,D);  Mtmin=Mtmax; 
mx=find(Tmax==-999);  mn=find(Tmin==-999); 
Mtmax(mx)=0;  Mtmin(mn)=0; 
Tmaxw(mx)=0;  Tmaxd(mx)=0;
Tminw(mn)=0;  Tmind(mn)=0;
%
% daily averages and standard deviations are calculated
%la fonction stats2.m est utilis?a la place de stats.m pour �viter les
%valeurs �gales a NaN.(Annie 2004)
v=ones(Y,1);
[aTmaxw,sTmaxw] = stats(Tmaxw,PnP.*Mtmax,v);
[aTminw,sTminw] = stats(Tminw,PnP.*Mtmin,v);
[aTmaxd,sTmaxd] = stats(Tmaxd,PnPm1.*Mtmax,v);
[aTmind,sTmind] = stats(Tmind,PnPm1.*Mtmin,v);
%
% smooth the series using Fourier transform with 2 harmonics (Richardson, 1981)
%
%  means
%
X=[ones(size(n)) sin(n/T)  cos(n/T) sin(2*n/T) cos(2*n/T)];
[aTmaxwest2,aC0(1),aC1(1),aC2(1),aD1(1),aD2(1)]=fourier(X,aTmaxw,n,T);
[aTminwest2,aC0(2),aC1(2),aC2(2),aD1(2),aD2(2)]=fourier(X,aTminw,n,T);
[aTmaxdest2,aC0(3),aC1(3),aC2(3),aD1(3),aD2(3)]=fourier(X,aTmaxd,n,T);
[aTmindest2,aC0(4),aC1(4),aC2(4),aD1(4),aD2(4)]=fourier(X,aTmind,n,T);
%
% standard deviations
%
[sTmaxwest2,sC0(1),sC1(1),sC2(1),sD1(1),sD2(1)]=fourier(X,sTmaxw,n,T);
[sTminwest2,sC0(2),sC1(2),sC2(2),sD1(2),sD2(2)]=fourier(X,sTminw,n,T);
[sTmaxdest2,sC0(3),sC1(3),sC2(3),sD1(3),sD2(3)]=fourier(X,sTmaxd,n,T);
[sTmindest2,sC0(4),sC1(4),sC2(4),sD1(4),sD2(4)]=fourier(X,sTmind,n,T);
%
%  plot the results
%
aTmaxwest2=aTmaxwest2';
aTmaxdest2=aTmaxdest2';
sTmaxwest2=sTmaxwest2';
sTmaxdest2=sTmaxdest2';
aTminwest2=aTminwest2';
aTmindest2=aTmindest2';
sTminwest2=sTminwest2';
sTmindest2=sTmindest2';

waitbar(0.98,h);
%
% calculate the residuals using the smoothed means and standard deviations
% one series of residuals per variables
%
[resTmax]=residual(Tmaxw,Tmaxd,aTmaxwest2,aTmaxdest2,sTmaxwest2,sTmaxdest2,Y,D);
[resTmin]=residual(Tminw,Tmind,aTminwest2,aTmindest2,sTminwest2,sTmindest2,Y,D);
%
% calculate the covariance matrix M0 (lag 0) and M1 (lag 1)
%
l=Y*D;
lag=0;
M0=covar2(resTmax,resTmin,lag,l);
lag=1;
M1=covar2(resTmax,resTmin,lag,l);
%
% calculate matrices A and C=BB' (Richardson (1981))
%
A=M1*inv(M0);
C=M0-M1*inv(M0)*M1';

% the correlation matrix has to be positive definite. When dealing with 
% measured data, observational errors, biases and missing data all 
% contaminate the true correlation and in some cases this results in 
% negative eigenvalues and the cholesky decomposition is not possible.
% thus, if the eigenvalues are not positive, simple set it close to zero.
[V,D] = eig(C);

i=find(D<0);      % find eigenvalues smaller than 0

if length(i)>=1 
        D(i)=0.000001;    % set the negative eigenvalues close to zero
        C=(V*D*inv(V));   % reconstruct matrix
end

%
%  produce a diagonal matrix lambda whose elements are eigenvalues of C
%  because C is symmetric (M0 is symmetric) and positive definite (elements
% in diagonal are positive), the eigenvalues will all be positive and real
%
% lambda=diag(eig(C));
[ee,dd]=eig(C); 
%
%  we have lambda = B'B (Matalas, 1967).  Therefore, B is obtained by matrix decomposition
%  can be used to obtain B
%
% B=chol(lambda);
B=ee*sqrt(dd)*ee';
waitbar(0.97,h);
%
%  results of the analysis are now saved for use by the generator

waitbar(1,h); 
close(h);  % close waitbar
