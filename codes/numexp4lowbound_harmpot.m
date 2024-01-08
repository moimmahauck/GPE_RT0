% Script to reproduce the numerical experiments of the paper 
% "Mixed finite elements for the Gross-Pitaevskii eigenvalue problem: 
% a priori error analysis and guaranteed lower energy bound" 
% by D. Gallistl, M. Hauck, Y. Liang, D. Peterseim

close all; clear all;

%% test of lower energy bound for harmonic potential
% parameters
L = 8;
kappa = 1000;
tol = 1e-8;
maxiter = 100;
tol_gs = 1e-4;
maxiter_gs = 20;

maxlevel = 7;
energy_gs = zeros(1,maxlevel);

% evaluate potential at coarse mesh
coarselevel = 3;
Tcoarse = refineMesh(getMesh('scaledSquareSymFKrot',L),coarselevel);
potential_eval_coarse = harmpotQ0(computeMids(Tcoarse),2^(coarselevel+1),L);

for level = 0:maxlevel-1
    fprintf('computing for h = 2^{-%d}\n',4+level);
    % mesh
    [T,P0] = refineMesh(Tcoarse,level);
    % assemble RT0 matrices
    [B,C,M] = assembleRT0(T);

    % prolongate coarse potential
    potential_eval = P0*potential_eval_coarse;
    
    % nonlinear eigenvalue solver
    % initial iteration
    mids = computeMids(T);
    u0 = (1-mids(:,1).^2/L^2).*(1-mids(:,2).^2/L^2);
    u0 = u0/sqrt(u0'*(M*u0));
    
    % [u,lambda,energy,normres] = Amethod(u0,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs);
    [u,lambda,energy,normres] = Jmethod(u0,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs);
    energy_gs(level+1) = energy(end);
end % for

%% plot the results 
% plot energies of the goundstate of different meshes
f1 = figure('position',[100,100,500,500]);
hs = sqrt(2)*(2.^-(4:10));
semilogx(hs,energy_gs); 
hold on;
semilogx(hs,energy_gs./(1+4*hs.^2.*pi^(-2).*energy_gs)); 
legend('$E_h$','$E_h^\mathrm{pp}$','interpreter','latex','FontSize', 14);
axis('square');
xticks(flip(hs));
xlim([hs(end),hs(1)]);
xticklabels({'2^{-10}','2^{-9}','2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}'});
xlabel('$h$','interpreter','latex','FontSize', 14);

exportgraphics(f1,'plots/lowbound_harmpot.png');

% plot energy difference in loglog plot
f2 = figure('position',[100,100,500,500]);
energy_ref = 6.054263127652939; % result of Q2-approximation with energy-adaptive riemannian optimization
hs = sqrt(2)*(2.^-(0:maxlevel-1));
loglog(hs,energy_ref-energy_gs); 
hold on; loglog(hs,energy_ref-energy_gs./(1+4*hs.^2.*pi^(-2).*energy_gs)); 
loglog(hs,hs.^2,'k:');
legend('$E^\mathrm{ref} - E_h$','$E^\mathrm{ref} - E_h^\mathrm{pp}$','slope 2','interpreter','latex','FontSize', 14);
axis('square'); 
xticks(flip(hs));
xlim([hs(end),hs(1)]);
xticklabels({'2^{-10}','2^{-9}','2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}'});
xlabel('$h$','interpreter','latex','FontSize', 14);

exportgraphics(f2,'plots/lowbound_randpotdifflog.png');

% plot potential
f3 = figure('position',[100,100,500,500]);
plotDiscFun(Tcoarse,potential_eval_coarse);
axis('square');
xticks(linspace(-L,L,5)); xlim([-L,L]);
yticks(linspace(-L,L,5)); ylim([-L,L]);
colorbar

exportgraphics(f3,'plots/plot_harmpot.png');

% plot solution
f4 = figure('position',[100,100,500,500]);
plotDiscFun(T,u,'edgecolor','none');
axis('square');
xticks(linspace(-L,L,5)); xlim([-L,L]);
yticks(linspace(-L,L,5)); ylim([-L,L]);
colorbar

exportgraphics(f4,'plots/groundstate_harmpot.png');

%% potential function
function y = harmpotQ0(x,N,L)
    N = round(N); 
    ax = repmat(1:N,N,1);
    idx4coords = [reshape(ax',N^2,1),reshape(ax,N^2,1)];
    h = 1/N;
    coords = h.*(idx4coords-1) + [h h]/2;
    coords = L*(2*coords-ones(1,2));
    val = .5*(coords(:,1).^2 + coords(:,2).^2);
    x = .5*(x + L*ones(1,2))/L;
    x(x == 0) = eps;
    y = val(sub2ind([N N],ceil(N*x(:,1)),ceil(N*x(:,2))));
end % for
