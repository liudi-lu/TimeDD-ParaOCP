clc;clear;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Problem setting ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 1e-1; gam = 10; Tend = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Discretization ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jx = 2^5; x = linspace(0, 1, Jx+1); hx = x(2) - x(1);
Jt = 2^5; t = linspace(0, Tend, Jt+1); ht = t(2) - t(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Solve sol exact ------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = BuildCNMatrix(hx, Jx, ht, Jt, nu, speye(Jx-1), sparse(Jx-1, Jx-1), ...
    gam*speye(Jx-1), speye(Jx-1));
F = BuildCNRHS(Jx, ht, Jt, yTarget(x(2:end-1), t'), y0(x(2:end-1)), ...
    gam*yTarget(x(2:end-1), t(end)));
U = M\F; U = reshape(U, Jx-1, 2*Jt+2);
y = U(:, 1:Jt+1); lam = U(:, Jt+2:end);

figure
mesh(t, x, [gx(t); U(:, 1:Jt+1); gx(t)]);
xlabel('$$t$$','interpreter','latex');
ylabel('$$x$$','interpreter','latex');
zlabel('$$y(x, t)$$','interpreter','latex');
set(gca, 'FontSize', 15); set(gca, 'linewidth', 1);

figure
mesh(t, x, [gx(t); U(:, Jt+2:end); gx(t)]/nu);
xlabel('$$t$$','interpreter','latex');
ylabel('$$x$$','interpreter','latex');
zlabel('$$u(x, t)$$','interpreter','latex');
set(gca, 'FontSize', 15); set(gca, 'linewidth', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------  Decompostion  ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
al = floor(2*(Jt+1)/5); Niter = 10;
errySD1 = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 1, 'SD1');
errySD1th = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 0.975, 'SD1');
errySD2 = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 1, 'SD2');
errySD3 = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 1, 'SD3');
errySD4 = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 1, 'SD4');
errySN1 = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 1, 'SN1');
errySN1th = SolveHeatOCP1dTimeAS(hx, Jx, x, ht, Jt, t, al, ...
    nu, gam, y, lam, Niter, 0.975, 'SN1');

figure
semilogy(1:Niter, errySD1/errySD1(1), '-*', ...
    1:Niter, errySD1th/errySD1th(1), '--+', ...
    1:Niter, errySD2/errySD2(1), '--h', ...
    1:Niter, errySD3/errySD3(1), '-d', ...
    1:Niter, errySD4/errySD4(1), '-.x', ...
    1:Niter, errySN1/errySN1(1), '--o', ...
    1:Niter, errySN1th/errySN1th(1), '-.s', ...
    'linewidth', 1.5, 'MarkerSize', 12);
xlim([1 Niter]); ylim([1e-15 1e5]);
xlabel('Iteration', 'interpreter', 'latex');
ylabel('Error', 'interpreter', 'latex');
legend({'SD$$_1$$', 'SD$$_1^{\theta}$$', ...
    'SD$$_2$$', 'SD$$_3$$', 'SD$$_4$$', ...
    'SN$$_1$$', 'SN$$_1^{\theta}$$'}, ...
    'interpreter', 'latex', 'location', 'best');
set(gca, 'FontSize', 20); set(gca, 'linewidth', 1.5);