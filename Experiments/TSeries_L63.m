function [r1,r2,r3,r4,C_uobs1,C_uobs2,C_uobs3,C_uobs4,C_vobs1,C_vobs2,C_vobs3,C_vobs4,delta_obs1,delta_obs2,delta_obs3,delta_obs4,P,Peps]=TSeries_L63();
s = 10; rh = 28; b = 8/3;                % the Lorenz system
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3) ...
          x(:,1).*x(:,2)-b*x(:,3)];
h = 0.01; n = 10; f = @(x) rk4(v,x,h,n);

epsilon=0.1;
veps = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3)+epsilon*x(:,1)...
          x(:,1).*x(:,2)-b*x(:,3)];
feps = @(x) rk4(veps,x,h,n);


c = [0.1 0.1 27]; r = [30 30 40];
t = Tree(c,r); depth=16;


xold=[1,1,1];
%variance=sqrt(h);
variance=0;
f_values=[];
tic
N=50000;
for i=1:N;
    if i<1000;
        xnew=f(xold);
        xold=xnew;
    else;
        xnew=f(xold);
        f_values=[f_values xnew'];
        t.insert(xnew', depth);
        xold=xnew;
    end
end
toc

Integration=f_values;

boxplot3(t); view(20,30); axis tight; axis square;
xlabel('x'); ylabel('y'); zlabel('z');

n = 50; x = linspace(-1,1,n)'; [XX,YY,ZZ] = meshgrid(x,x,x);
X = [ XX(:) YY(:) ZZ(:) ];
X = 2*rand(100,3)-1;      
P = tpmatrix_TS(t, f,Integration', X, depth,1);
[v,lambda] = eigs(P,51,'lm');             % compute eigenvector at eigenvalue 1
vp = log10(max(abs(v(:,1)),1e-16)) ;
usub1=v;


plot(log(diag(lambda)),'.','markersize',10)

xold=[1,1,1];
feps_values=[];
for i=1:N;
    if i<1000;
        xnew=feps(xold);
        xold=xnew;
    else;
        xnew=feps(xold);
        feps_values=[feps_values xnew'];
        xold=xnew;
    end
end


n = 50; x = linspace(-1,1,n)'; [XX,YY,ZZ] = meshgrid(x,x,x);
X = [ XX(:) YY(:) ZZ(:) ];
X = 2*rand(100,3)-1;      
Peps = tpmatrix_TS(t, feps,Integration', X, depth,1);
[v,lambda] = eigs(Peps,1,'lm');             % compute eigenvector at eigenvalue 1
vp = log10(max(abs(v(:,1)),1e-16)) ;
vsub1=v;


b=t.boxes(-1);
d = t.dim;
N = size(b,2);
K=[1:N];
c = b(1:d,K);

%Observables
obs1=c(1,:).^2;
obs2=c(2,:).^2;
obs3=c(3,:).^2;
obs4=c(3,:);


C_uobs1=norm(obs1'.*usub1,1)/(norm(usub1,1));
C_uobs2=norm(obs2'.*usub1,1)/(norm(usub1,1));
C_uobs3=norm(obs3'.*usub1,1)/(norm(usub1,1));
C_uobs4=norm(obs4'.*usub1,1)/(norm(usub1,1));


C_vobs1=norm(obs1'.*vsub1,1)/(norm(vsub1,1));
C_vobs2=norm(obs2'.*vsub1,1)/(norm(vsub1,1));
C_vobs3=norm(obs3'.*vsub1,1)/(norm(vsub1,1));
C_vobs4=norm(obs4'.*vsub1,1)/(norm(vsub1,1));

m=(Peps-P);
%[m,mmm]=lin_resp(P);
tic
psi1=Response_Operator(P,m);
toc

psi1obs1=psi1*obs1';
psi1obs2=psi1*obs2';
psi1obs3=psi1*obs3';
psi1obs4=psi1*obs4';

delta_obs1=norm(usub1.*psi1obs1,1)/(norm(usub1,1));
delta_obs2=norm(usub1.*psi1obs2,1)/(norm(usub1,1));
delta_obs3=norm(usub1.*psi1obs3,1)/(norm(usub1,1));
delta_obs4=norm(usub1.*psi1obs4,1)/(norm(usub1,1));

r1=10*abs(C_uobs1 - C_vobs1);
r2=10*abs(C_uobs2 - C_vobs2);
r3=10*abs(C_uobs3 - C_vobs3);
r4=10*abs(C_uobs4 - C_vobs4);
    
% x=Integration(1,:);
% N = length(x);
% Fs=N;
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% f2 = 0:Fs/length(x):Fs/2;
% semilogx(f2,10*log10(psdx)); 
% plot(f2,psdx)
% 
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

end

