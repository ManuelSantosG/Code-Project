tic
[M,usub1,t,depth]=Complete(0,12);
[Meps,vsub1,t,depth]=Complete(0.1,12);
m=Meps-M;
toc

tic
psi1=Response_Operator(M,m);
toc

b=t.boxes(-1);
d = t.dim;
N = size(b,2);
K=[1:N];
c = b(1:d,K); %Centers of the boxes

%Observables
obs1=c(1,:).^2;
obs2=c(2,:).^2;
obs3=c(3,:).^2;
obs4=c(3,:);

psi1obs1=psi1*obs1';
psi1obs2=psi1*obs2';
psi1obs3=psi1*obs3';
psi1obs4=psi1*obs4';


C_uobs1=norm(obs1'.*usub1,1)/(norm(usub1,1));
C_uobs2=norm(obs2'.*usub1,1)/(norm(usub1,1));
C_uobs3=norm(obs3'.*usub1,1)/(norm(usub1,1));
C_uobs4=norm(obs4'.*usub1,1)/(norm(usub1,1));

delta_obs1=norm(psi1obs1'.*usub1,1)*0.1/(norm(usub1,1));
delta_obs2=norm(psi1obs2'.*usub1,1)*0.1/(norm(usub1,1));
delta_obs3=norm(psi1obs3'.*usub1,1)*0.1/(norm(usub1,1));
delta_obs4=norm(psi1obs4'.*usub1,1)*0.1/(norm(usub1,1));

C_vobs1=norm(obs1'.*vsub1,1)/(norm(vsub1,1));
C_vobs2=norm(obs2'.*vsub1,1)/(norm(vsub1,1));
C_vobs3=norm(obs3'.*vsub1,1)/(norm(vsub1,1));
C_vobs4=norm(obs4'.*vsub1,1)/(norm(vsub1,1));

r1=10*abs(C_uobs1 - C_vobs1)
r2=abs(C_uobs2 - C_vobs2)
r3=10*abs(C_uobs3 - C_vobs3)
r4=10*abs(C_uobs4 - C_vobs4)
