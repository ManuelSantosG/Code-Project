


B = chebop(-1, 1);
B.op = @(x,u)   diff(u,2) - x.*diff(u,1) - u;
B.lbc=0;
B.rbc=0;

lam = eigs(B,100);
MS = 'markersize';
clf, plot(lam, 'r.', MS, 16), grid on, axis equal
spectral_abscissa = max(real(lam))