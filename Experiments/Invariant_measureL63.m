% GAIO demo: Natural invariant measure of the logistic map

s = 10; rh = 28; b = 0.4;                % the Lorenz system
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3) ...
          x(:,1).*x(:,2)-b*x(:,3)];

h = 0.01; n = 10; f = @(x) rk4(v,x,h,n); % f is the flow     
      
      
n = 300; X1 = linspace(-1,1,n)'; X = [X1,X1,X1];        % uniform grid of sample points
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

% Plot invariant density

figure
vp = log10(max(abs(v(:,1)),1e-16)) ;
boxplot3(t,'depth',depth,'density', vp ,'alpha',0.1);

% n = t.count(depth); x = linspace(0,1,n);
% h = abs(v)*n./norm(v,1);
% bar(x,h,1); axis([0 1 0 10]);
% xlabel('x'); ylabel('density');


% Published with MATLAB® R2013a