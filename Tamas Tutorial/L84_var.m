function f = L84_var(t,X,a,b,F,G,Gp)

% 30.10.17. RHS of ODE and variational eq. coupled for L84.

% State variables
x = X(1); y = X(2); z = X(3);

% Components of the Jacobian
Y = [X(4), X(7), X(10);
     X(5), X(8), X(11);
     X(6), X(9), X(12)];

f=zeros(12,1);

% RHS of L84 
f(1 ) = -(y^2 + z^2)  - a*x + a*F;
f(2 ) =   x*y - b*x*z -   y + G  ;
f(3 ) =   x*z + b*x*y -   z + Gp ;

% Jacobian
J = zeros(3);

% Use linear indexing
J(1) = -a;
J(2) = y - b*z;
J(3) = b*y + z;

J(4) = -2*y;
J(5) = x - 1;
J(6) = b*x;

J(7) = -2*z;
J(8) = -b*x;
J(9) = x - 1;

% Variational equation   
f(4:12)=J*Y;