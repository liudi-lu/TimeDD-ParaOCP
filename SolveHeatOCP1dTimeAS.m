function [gy, glam] = SolveHeatOCP1dTimeAS(h, J, x, dt, N, t, al, ...
    nu, gam, y, lam, Niter, th, Var)
% Alternating Schwarz to solve the HeatOCP with two subdomains I1 and I2.
% We construct the matrix Mj globally and the right hand side Fj at 
% each iteration, since Fj depends on the transmission conditions 
% gl and gr.
%
% Parameters
% ----------
% h : float
%     space meshsize
% J : int
%     number of points in space (J+1)
% x : vector
%     grid in space
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% t : vector
%     grid in time
% al : int
%     index of the interface alpha in time
% nu : float
%     penalization for control
% gam : float
%     penalization for final target
% y : matrix
%     exact solution of state variable y
% lam : matrix
%     exact solution of adjoint variable lambda
% Niter : int
%     number of iteration
% th : float
%     relaxation parameter
% Var : string
%     name of variant  

% Initial guess of the transmission condition
gl = zeros(J-1, 1); 

% Error of y and lambda
gy = zeros(1, Niter); glam = gy;

% construct matrices of each subdomain
if strcmp(Var, 'SD1')
    M1 = BuildCNMatrix(h, J, dt, al, nu, speye(J-1), ...
        sparse(J-1, J-1), sparse(J-1, J-1), speye(J-1));
    M2 = BuildCNMatrix(h, J, dt, N-al, nu, speye(J-1), ...
        sparse(J-1, J-1), gam*speye(J-1), speye(J-1));
elseif strcmp(Var, 'SD2')
    M1 = BuildCNMatrix(h, J, dt, al, nu, speye(J-1), ...
        sparse(J-1, J-1), speye(J-1), sparse(J-1, J-1));
    M2 = BuildCNMatrix(h, J, dt, N-al, nu, sparse(J-1, J-1), ...
        speye(J-1), gam*speye(J-1), speye(J-1));
elseif strcmp(Var, 'SD3')
    M1 = BuildCNMatrix(h, J, dt, al, nu, speye(J-1), ...
        sparse(J-1, J-1), speye(J-1), sparse(J-1, J-1));
    M2 = BuildCNMatrix(h, J, dt, N-al, nu, speye(J-1), ...
        sparse(J-1, J-1), gam*speye(J-1), speye(J-1));
elseif strcmp(Var, 'SD4')
    M1 = BuildCNMatrix(h, J, dt, al, nu, speye(J-1), ...
        sparse(J-1, J-1), sparse(J-1, J-1), speye(J-1));
    M2 = BuildCNMatrix(h, J, dt, N-al, nu, sparse(J-1, J-1), ...
        speye(J-1), gam*speye(J-1), speye(J-1));
elseif strcmp(Var, 'SN1')
    Ax = A1d(0, h, J-1);
    M1 = BuildCNMatrix(h, J, dt, al, nu, speye(J-1), ...
        sparse(J-1, J-1), speye(J-1), Ax);
    M2 = BuildCNMatrix(h, J, dt, N-al, nu, -Ax, speye(J-1)/nu, ...
        gam*speye(J-1), speye(J-1));
else
    disp('Error in the name of variant !!!')
    return
end

% Alternating Schwarz variants' iterations
for n = 1 : Niter
    %construct rhs in I1
    F1 = BuildCNRHS(J, dt, al, yTarget(x(2:end-1), t(1:al+1)'), ...
        y0(x(2:end-1)), gl);
    % solve the linear system in I1
    U1 = M1\F1; U1 = reshape(U1, J-1, 2*al+2);
    y1 = U1(:, 1:al+1); lam1 = U1(:, al+2:end);
    
    % pass transmission condition to I2
    if strcmp(Var, 'SD1') || strcmp(Var, 'SD3') 
        gr = y1(:, end);
    elseif strcmp(Var, 'SD2') || strcmp(Var, 'SD4')
        gr = lam1(:, end);
    elseif strcmp(Var, 'SN1')
        gr = -Ax*y1(:, end) + lam1(:, end)/nu;
    else
        disp('Error in the name of variant !!!')
        return
    end

    %construct rhs in I2
    F2 = BuildCNRHS(J, dt, N-al, yTarget(x(2:end-1), t(al+1:end)'), ...
        gr, gam*yTarget(x(2:end-1), t(end)));
    %solve the linear system in I2
    U2 = M2\F2; U2 = reshape(U2, J-1, 2*(N-al)+2); 
    y2 = U2(:, 1:N-al+1); lam2 = U2(:, N-al+2:end);

    %update transmission condition in I1
    if strcmp(Var, 'SD1') || strcmp(Var, 'SD4') 
        gl = (1-th)*gl + th*lam2(:, 1);
    elseif strcmp(Var, 'SD2') || strcmp(Var, 'SD3')
        gl = (1-th)*gl + th*y2(:, 1);
    elseif strcmp(Var, 'SN1')
        gl = (1-th)*gl + th*(y2(:, 1) + Ax*lam2(:, 1));
    else
        disp('Error in the name of variant !!!')
        return
    end

    %compute error with infinity norm
    gy(n) = norm(y - [y1(:, 1:end-1) y2], 'inf');
    glam(n) = norm(lam - [lam1 lam2(:, 2:end)], 'inf');
end
end