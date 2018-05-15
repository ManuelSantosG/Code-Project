% GAIO demo: Natural invariant measure of the logistic map
function [M,v,t,depth]=Complete(epsilon,depth);
    s = 10; rh = 28; b = 8/3;                % the Lorenz system
    v = @(x) [s*(x(:,2)-x(:,1)) ...
              rh*x(:,1)-x(:,2)-x(:,1).*x(:,3)...
              x(:,1).*x(:,2)-b*x(:,3)];
    veps = @(x) [s*(x(:,2)-x(:,1)) ...
          rh*x(:,1)-x(:,2)-x(:,1).*x(:,3)+epsilon*x(:,1)...
          x(:,1).*x(:,2)-b*x(:,3)];
    h = 0.01; n = 10; f = @(x) rk4(v,x,h,n); % f is the flow
    feps= @(x) rk4(veps,x,h,n);

    n = 100; x = linspace(-1,1,n)'; [XX,YY,ZZ] = meshgrid(x,x,x);
    X = [ XX(:) YY(:) ZZ(:) ];               % $7^3$ sample points
    c = [0.1 0.1 27]; r = [30 30 40];
    t = Tree(c,r);                           % the box collection


    a = sqrt(b*(rh-1)); x0 = [a a rh-1;-a -a rh-1]; % equilibria 
    t.insert(x0', depth);         % construct [x0]
    gum(t, f, X, depth);                     % compute global unstable manifold
    
    boxplot3(t); view(20,30); axis tight; axis square;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    X = 2*rand(100,3)-1;                     % points for Monte Carlo quadrature
    P = tpmatrix(t, feps, X, depth,1);% transition matrix
    [v,lambda] = eigs(P,1,'lm');             % compute eigenvector at eigenvalue 1
    vp = log10(max(abs(v(:,1)),1e-16)) ;     % plot measure logarithmically
    
%     figure
%     boxplot3(t,'depth',depth,'density', vp ,'alpha',0.1);
    
    %M=full(P);
    M=P;
    
    
    %Plot invariant density

    figure
    vp = log10(max(abs(v(:,1)),1e-16)) ;
    boxplot3(t,'depth',depth,'density', vp ,'alpha',0.1);

    % n = t.count(depth); x = linspace(0,1,n);
    % h = abs(v)*n./norm(v,1);
    % bar(x,h,1); axis([0 1 0 10]);
    % xlabel('x'); ylabel('density');


    % Published with MATLAB® R2013a
end