function M = BuildCNMatrix(h, J, dt, N, nu, Ay0, Alam0, AyT, AlamT)
% Build the matrix M to solve the forward-backward linear system  
% MU = F, with U = [Y; Lambda]. There are four components in the matrix 
% M = [M11, M12; M21, M22]. The first block of M11 and M12 and the last
% block of M21 and M22 are specially adapted to different transmission
% conditions. The numerical scheme is Crank-Nicolson.
%
% Parameters
% ----------
% h : float
%     space meshsize
% J : int
%     number of points in space (J+1)
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% nu : float
%     penalization for control
% Ay0 : matrix
%     a (J-1)*(J-1) matrix for the first block of M11
% Alam0 : matrix
%     a (J-1)*(J-1) matrix for the first block of M12
% AyT : matrix
%     a (J-1)*(J-1) matrix for the last block of M21
% AlamT : matrix
%     a (J-1)*(J-1) matrix for the last block of M22

% generate the spatial matrix
A = A1d(0, h, J-1);
e = ones(N+1, 1);

% construct four matrices
M11 = kron(speye(N+1), speye(J-1)+dt*A/2) - ...
    kron(spdiags(e, -1, N+1, N+1), speye(J-1)-dt*A/2);
M22 = -M11';
M21 = kron(spdiags([e e], [0 1], N+1, N+1), -dt*speye(J-1)/2);
M12 = M21'/nu;

% adapt the initial and final condition
M11(1:J-1, 1:J-1) = Ay0;
M12(1:J-1, 1:J-1) = Alam0;
M21(end-J+2:end, end-J+2:end) = AyT;
M22(end-J+2:end, end-J+2:end) = AlamT;

% concatenation of four matrices
M = [M11, M12; M21, M22];
end