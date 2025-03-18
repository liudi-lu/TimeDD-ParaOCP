function g = yTarget(x, t)
% Target function \hat y(x, t), which satisfies the initial condition
% \hat y(x, 0) = 0, and two Dirichlet boundary conditions 
% \hat y(0, t) = 0 and \hat y(1, t) = 0 for x in [0, 1].
g = sin(pi*x).*(2*t.^2+t);
end