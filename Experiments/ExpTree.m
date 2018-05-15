s = 10; rh = 28; b = 0.4;                % the Lorenz system
v = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3) ...
          x(:,1).*x(:,2)-b*x(:,3)];
h = 0.01; n = 10; f = @(x) rk4(v,x,h,n); % f is the flow


n = 7; x = linspace(-1,1,n)'; [XX,YY,ZZ] = meshgrid(x,x,x);
X = [ XX(:) YY(:) ZZ(:) ];               % $7^3$ sample points
c = [0.1 0.1 27]; r = [30 30 40];
t = Tree(c,r);                           % the box collection


a = sqrt(b*(rh-1)); x0 = [a a rh-1;-a -a rh-1]; % equilibria
depth = 15; 
t.insert(x0', depth);         % construct [x0]
gum(t, f, X, depth);                     % compute global unstable manifold

boxplot3(t); view(20,30); axis tight; axis square;
xlabel('x'); ylabel('y'); zlabel('z');

