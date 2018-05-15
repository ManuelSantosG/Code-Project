% GAIO demo: Natural invariant measure of the logistic map

s = 10; rh = 28; b = 8/3;                % the Lorenz system
      
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3) ...
          x(:,1).*x(:,2)-b*x(:,3)];

h = 0.01; n = 10; f = @(x) rk4(v,x,h,n); % f is the flow     
      
n = 100; X1 = linspace(-1,1,n)'; X = [X1,X1,X1];        % uniform grid of sample points
c = [0.1,0.1,27]; r = [30,30,40]; t = Tree(c, r);  % the tree
% Construct full subdivison
sd = 8; depth = 15;
for i=1:depth
    t.set_flags('all', sd);
    t.subdivide;
end
% Compute invariant vector




P = tpmatrix(t, f, X, depth,1);
[v,lambda] = eigs(P,1);
%PP=full(P);


b=t.boxes(depth);
d = t.dim;
N = size(b,2);
K=[1:N];
c = b(1:d,K);


Ic=t.search(c,depth);
%Ic=sort(Ic);

fff=c(1,:).^2;
fff1=c(2,:).^2;
fff2=c(3,:).^2;
fff3=c(3,:);


NN=10000;

%C=eye(size(PP))/NN;
%size(C)
%PPi_1=eye(size(PP));
% %C=fff';
% for i =1:NN-1
%     PPi=PP*PPi_1;
%     C=C+PPi/NN;
%     PPi_1=PPi;
% %     plot(C)
% %     hold on
% end
% %C1=norm(C,1);
% C1=C*fff';
C_xsq=norm(fff'.*v,1)/(norm(v,1));
C_ysq=norm(fff1'.*v,1)/(norm(v,1));
C_zsq=norm(fff2'.*v,1)/(norm(v,1));
C_z=norm(fff3'.*v,1)/(norm(v,1));
% Plot invariant density

% figure
% vp = log10(max(abs(v(:,1)),1e-16)) ;
% boxplot3(t,'depth',depth,'density', vp ,'alpha',0.1);

% n = t.count(depth); x = linspace(0,1,n);
% h = abs(v)*n./norm(v,1);
% bar(x,h,1); axis([0 1 0 10]);
% xlabel('x'); ylabel('density');


% Published with MATLAB® R2013a