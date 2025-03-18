function A = A1d(eta, h, J)
% Finite Difference approximation of the one dimensional operator 
% (-Delta + eta) for the meshsize h and J points.
e = ones(J, 1);
A = spdiags([-e/h^2 (eta+2/h^2)*e -e/h^2], [-1 0 1], J, J);
end