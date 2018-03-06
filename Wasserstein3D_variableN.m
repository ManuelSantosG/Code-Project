clear all
close all
   
getd = @(p)path(p,path);       %Computation of optimal transport cost - needed for Gamma matrix
getd('toolbox_signal/');
getd('toolbox_general/');
flat = @(x)x(:);
Cols = @(n0,n1)sparse( flat(repmat(1:n1, [n0 1])), ...
             flat(reshape(1:n0*n1,n0,n1) ), ...
             ones(n0*n1,1) );
Rows = @(n0,n1)sparse( flat(repmat(1:n0, [n1 1])), ...
             flat(reshape(1:n0*n1,n0,n1)' ), ...
             ones(n0*n1,1) );
Sigma = @(n0,n1)[Rows(n0,n1);Cols(n0,n1)];

maxit = 1e20; tol = 1e-9;
otransp = @(C,p0,p1)reshape( perform_linprog( ...
        Sigma(length(p0),length(p1)), ...
        [p0(:);p1(:)], C(:), 0, maxit, tol), [length(p0) length(p1)] ); 
        
minX=-2.0;    % Here I set the range for data binning
rangeX=6.0;
minY=-3.5;
rangeY=7.0;
minZ=-3.5;
rangeZ=7.0;
% minX=-2.5;
% rangeX=7.0;
% minY=-4.0;
% rangeY=8.0;
% minZ=-4.0;
% rangeZ=8.0;

load 'L84stat.d';     % Loading files
load 'L84_63stat.d';
load 'L84_1stat.d';
load 'L84_2stat.d';

min(L84_63stat(:,1))     % This is just to check the extremes of the n dimensional cube in which apply data binning
min(L84_2stat(:,1))
max(L84_63stat(:,1))
max(L84_2stat(:,1))
min(L84_63stat(:,2))
min(L84_2stat(:,2))
max(L84_63stat(:,2))
max(L84_2stat(:,2))
min(L84_63stat(:,3))
min(L84_2stat(:,3))
max(L84_63stat(:,3))
max(L84_2stat(:,3))

fileB = fopen('exp3Dbasic.txt','a');     %Open files where I write the Wasserstein distances
file1 = fopen('exp3Dorder1.txt','a');
file2 = fopen('exp3Dorder2.txt','a');

format shortg

for numdiv=1:10   % Loop on number of cubes per side
  
sample0=L84_63stat';   
sample1=L84stat';

totlen=numdiv^3;     %Number of possible coordinates
supp0=zeros(numdiv,numdiv,numdiv);   %Array 3D where I put how many points there are for each little cube
supp1=zeros(numdiv,numdiv,numdiv);
X0=zeros(3,totlen);    %Array with possible coordinates
X1=zeros(3,totlen);
p0=zeros(totlen,1);    %Array with weights of single coordinates
p1=zeros(totlen,1);
for c=1:length(sample0)
  for i=1:numdiv
    for j=1:numdiv
      for k=1:numdiv
        if (sample0(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample0(1,c)< (minX+rangeX*i/numdiv))
          if (sample0(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample0(2,c)< (minY+rangeY*j/numdiv))
            if (sample0(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample0(3,c)< (minZ+rangeZ*k/numdiv))
              supp0(i,j,k)=supp0(i,j,k)+1;     %Increasing number of points per cube
            end
          end
        end
        if (sample1(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample1(1,c)< (minX+rangeX*i/numdiv))
          if (sample1(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample1(2,c)< (minY+rangeY*j/numdiv))
            if (sample1(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample1(3,c)< (minZ+rangeZ*k/numdiv))
              supp1(i,j,k)=supp1(i,j,k)+1;
            end
          end
        end
      end
    end
  end
end  
c=0;
for i=1:numdiv
  for j=1:numdiv
     for k=1:numdiv
        c=c+1;
        X0(1,c)=i;   %Array with coordinates of cubes
        X0(2,c)=j;
        X0(3,c)=k;
     end
  end
end
X1=X0;
normalize = @(a)a/sum(a(:));
supp0 = normalize(supp0);    %Normalization
supp1 = normalize(supp1);
for c=1:totlen
  p0(c)=supp0(X0(1,c),X0(2,c),X0(3,c));    %Density of points array
  p1(c)=supp1(X1(1,c),X1(2,c),X1(3,c));
end

% Deleting the empty cubes and their null density - it speeds up
% computations
arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p0(c)~=0)
        arrayuno=cat(1,arrayuno,p0(c));
        arraydue=cat(2,arraydue,X0(:,c));
    end
end

p0=arrayuno;
X0=arraydue;

arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p1(c)~=0)
        arrayuno=cat(1,arrayuno,p1(c));
        arraydue=cat(2,arraydue,X1(:,c));
    end
end

p1=arrayuno;
X1=arraydue;

C = repmat( sum(X0.^2)', [1 length(p1)] ) + ...      %C matrix
    repmat( sum(X1.^2), [length(p0) 1] ) - 2*X0'*X1;

gamma = otransp(C,p0,p1);              %Optimal Gamma matrix

wasserstein=0.;
for i=1:length(p0)
    for j=1:length(p1)
        wasserstein=wasserstein+gamma(i,j)*C(i,j);    %Wasserstein^2 computation
    end
end
wasserstein_basic=sqrt(wasserstein);              %Wasserstein distance


sample0=L84_63stat';
sample1=L84_1stat';

totlen=numdiv^3;
supp0=zeros(numdiv,numdiv,numdiv);
supp1=zeros(numdiv,numdiv,numdiv);
X0=zeros(3,totlen);
X1=zeros(3,totlen);
p0=zeros(totlen,1);
p1=zeros(totlen,1);
for c=1:length(sample0)
  for i=1:numdiv
    for j=1:numdiv
      for k=1:numdiv
        if (sample0(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample0(1,c)< (minX+rangeX*i/numdiv))
          if (sample0(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample0(2,c)< (minY+rangeY*j/numdiv))
            if (sample0(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample0(3,c)< (minZ+rangeZ*k/numdiv))
              supp0(i,j,k)=supp0(i,j,k)+1;
            end
          end
        end
        if (sample1(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample1(1,c)< (minX+rangeX*i/numdiv))
          if (sample1(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample1(2,c)< (minY+rangeY*j/numdiv))
            if (sample1(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample1(3,c)< (minZ+rangeZ*k/numdiv))
              supp1(i,j,k)=supp1(i,j,k)+1;
            end
          end
        end
      end
    end
  end
end  
c=0;
for i=1:numdiv
  for j=1:numdiv
     for k=1:numdiv
        c=c+1;
        X0(1,c)=i; 
        X0(2,c)=j;
        X0(3,c)=k;
     end
  end
end
X1=X0;
normalize = @(a)a/sum(a(:));
supp0 = normalize(supp0); 
supp1 = normalize(supp1);
for c=1:totlen
  p0(c)=supp0(X0(1,c),X0(2,c),X0(3,c)); 
  p1(c)=supp1(X1(1,c),X1(2,c),X1(3,c));
end

arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p0(c)~=0)
        arrayuno=cat(1,arrayuno,p0(c));
        arraydue=cat(2,arraydue,X0(:,c));
    end
end

p0=arrayuno;
X0=arraydue;

arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p1(c)~=0)
        arrayuno=cat(1,arrayuno,p1(c));
        arraydue=cat(2,arraydue,X1(:,c));
    end
end

p1=arrayuno;
X1=arraydue;

C = repmat( sum(X0.^2)', [1 length(p1)] ) + ...
    repmat( sum(X1.^2), [length(p0) 1] ) - 2*X0'*X1;

gamma = otransp(C,p0,p1);  

wasserstein=0.;
for i=1:length(p0)
    for j=1:length(p1)
        wasserstein=wasserstein+gamma(i,j)*C(i,j); 
    end
end
wasserstein_1=sqrt(wasserstein); 


sample0=L84_63stat';
sample1=L84_2stat';

totlen=numdiv^3;  
supp0=zeros(numdiv,numdiv,numdiv); 
supp1=zeros(numdiv,numdiv,numdiv);
X0=zeros(3,totlen); 
X1=zeros(3,totlen);
p0=zeros(totlen,1);  
p1=zeros(totlen,1);
for c=1:length(sample0)
  for i=1:numdiv
    for j=1:numdiv
      for k=1:numdiv
        if (sample0(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample0(1,c)< (minX+rangeX*i/numdiv))
          if (sample0(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample0(2,c)< (minY+rangeY*j/numdiv))
            if (sample0(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample0(3,c)< (minZ+rangeZ*k/numdiv))
              supp0(i,j,k)=supp0(i,j,k)+1; 
            end
          end
        end
        if (sample1(1,c)>=(minX+rangeX*(i-1)/numdiv)) && (sample1(1,c)< (minX+rangeX*i/numdiv))
          if (sample1(2,c)>=(minY+rangeY*(j-1)/numdiv)) && (sample1(2,c)< (minY+rangeY*j/numdiv))
            if (sample1(3,c)>=(minZ+rangeZ*(k-1)/numdiv)) && (sample1(3,c)< (minZ+rangeZ*k/numdiv))
              supp1(i,j,k)=supp1(i,j,k)+1;
            end
          end
        end
      end
    end
  end
end  
c=0;
for i=1:numdiv
  for j=1:numdiv
     for k=1:numdiv
        c=c+1;
        X0(1,c)=i; 
        X0(2,c)=j;
        X0(3,c)=k;
     end
  end
end
X1=X0;
normalize = @(a)a/sum(a(:));
supp0 = normalize(supp0);   
supp1 = normalize(supp1);
for c=1:totlen
  p0(c)=supp0(X0(1,c),X0(2,c),X0(3,c));  
  p1(c)=supp1(X1(1,c),X1(2,c),X1(3,c));
end

arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p0(c)~=0)
        arrayuno=cat(1,arrayuno,p0(c));
        arraydue=cat(2,arraydue,X0(:,c));
    end
end

p0=arrayuno;
X0=arraydue;

arrayuno=[];
arraydue=[];
for c=1:totlen
    if (p1(c)~=0)
        arrayuno=cat(1,arrayuno,p1(c));
        arraydue=cat(2,arraydue,X1(:,c));
    end
end

p1=arrayuno;
X1=arraydue;

C = repmat( sum(X0.^2)', [1 length(p1)] ) + ... 
    repmat( sum(X1.^2), [length(p0) 1] ) - 2*X0'*X1;

gamma = otransp(C,p0,p1);    

wasserstein=0.;
for i=1:length(p0)
    for j=1:length(p1)
        wasserstein=wasserstein+gamma(i,j)*C(i,j);  
    end
end
wasserstein_2=sqrt(wasserstein);   

% Printing on screen and on file the results - WD is rescaled
clock
fprintf('Number of cubes per side= %d\n', numdiv);
fprintf('Total number of cubes= %d\n', numdiv^3);
fprintf('Wasserstein distance - basic: %f\n', wasserstein_basic*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);
fprintf('Wasserstein distance - order1: %f\n', wasserstein_1*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);
fprintf('Wasserstein distance - order2: %f\n\n', wasserstein_2*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);

fprintf(fileB,'%f   ', wasserstein_basic*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);
fprintf(file1,'%f   ', wasserstein_1*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);
fprintf(file2,'%f   ', wasserstein_2*((rangeX*rangeY*rangeZ)^(1./3.))/numdiv);

end
fprintf(fileB,'\n');
fprintf(file1,'\n');
fprintf(file2,'\n');
fclose(fileB);
fclose(file1);
fclose(file2);