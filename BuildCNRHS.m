function F = BuildCNRHS(J, dt, N, Yhat, F0, FT)
% Build the vector F to solve the problem MU = F. There are two 
% components in the vector F = [F1; F2], where F1 corresponds to 
% the right hand side of the state variable Y, and F2 corresponds 
% to the right hand side of the adjoint variable Lambda.
%
% Parameters
% ----------
% J : int
%     number of points in space (J+1)
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% Yhat : matrix
%     target function \hat y(x, t)
% F0 : vector
%     initial condition
% FT : vector
%     final condition

% construct F1 to solve Y
F1 = zeros((N+1)*(J-1), 1);

% add the initial condition
F1(1:J-1) = F0;

% construct F2 to solve Lambda
F2 = -dt*(Yhat(2:end, :) + Yhat(1:end-1, :))'/2;

% add the final condition
F2 = [F2(:); FT(:)];

% concatenate F1 and F2
F = [F1; F2];
end