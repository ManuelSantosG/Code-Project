function yn = RK4(odefun,ra,h,ya,varargin) 
% Classical 4th order Runge-Kutta solver with step size h -- 
% takes only a single step


neq = length(ya);
F = zeros(neq,4);

F(:,1) = feval(odefun,ra,ya,varargin{:});
F(:,2) = feval(odefun,ra + h / 2,ya + h / 2 * F(:,1),varargin{:});
F(:,3) = feval(odefun,ra + h / 2,ya + h / 2 * F(:,2),varargin{:});  
F(:,4) = feval(odefun,ra + h,ya + h * F(:,3),varargin{:});

yn = ya + (h / 6) * (F(:,1) + 2 * F(:,2) + 2 * F(:,3) + F(:,4)); 