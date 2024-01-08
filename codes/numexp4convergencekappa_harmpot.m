% Script to reproduce the numerical experiments of the paper 
% "Mixed finite elements for the Gross-Pitaevskii eigenvalue problem: 
% a priori error analysis and guaranteed lower energy bound" 
% by D. Gallistl, M. Hauck, Y. Liang, D. Peterseim

close all; clear all;

%% convergence test for harmonic potential
% parameters
L = 8;
kappas = [1 10 100 1000];
potential = @(x) .5*(x(:,1).^2 + x(:,2).^2);
tol = 1e-8;
maxiter = 100;
tol_gs = 1e-4;
maxiter_gs = 20;

maxlevel = 7;
nref = 9;

%% compute reference solution on level nref
% mesh
Tref = refineMesh(getMesh('scaledSquareSymFKrot',L),nref);
% assemble RT0 matrices
[edgesref,nodes2edgesref] = getEdgeProperties(Tref);
[Bref,Cref,Mref] = assembleRT0(Tref);

% evaluate coarse potential
potential_evalref = potential(computeMids(Tref));

% nonlinear eigenvalue solver
% initial iteration
mids = computeMids(Tref);
u0ref = (1-mids(:,1).^2/L^2).*(1-mids(:,2).^2/L^2);
u0ref = u0ref/sqrt(u0ref'*(Mref*u0ref));

u_err = zeros(length(kappas),maxlevel+1);
gu_err = zeros(length(kappas),maxlevel+1);
lambda_err = zeros(length(kappas),maxlevel+1);
energy_err = zeros(length(kappas),maxlevel+1);

for indkappa = 1:length(kappas)
    kappa = kappas(indkappa);
    fprintf('computing for kappa = %d\n',kappa);
    fprintf('computing reference solution\n');

    % [u,lambda,energy,normres] = Amethod(u0,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs);
    [uref,lambdaref,energyref,normresref] = Jmethod(u0ref,Bref,Cref,Mref,potential_evalref,kappa,tol,maxiter,tol_gs,maxiter_gs);
    [guref,~,~,~] = pcg(Bref,-Cref'*uref);
    
    % reconstruct flux from vector guref
    guref_x = zeros(Tref.nelems,3);
    guref_y = zeros(Tref.nelems,3);
    for elem = 1:Tref.nelems
        % evaluate RT0 basis functions at local nodes
        globnodes = Tref.elems(elem,:);
        globedgeinds = full(nodes2edgesref(sub2ind([Tref.nnodes,Tref.nnodes],globnodes([2 1 1]),globnodes([3 3 2]))));
        loccoords = Tref.coords(globnodes,:)';
        voledges = [norm(loccoords(:,2)-loccoords(:,3),2), ...
                    norm(loccoords(:,1)-loccoords(:,3),2),...
                    norm(loccoords(:,1)-loccoords(:,2),2)];
        signs = -1 + 2*(edgesref(globedgeinds,3) == elem)';
        voledgeswsigns = diag(voledges.*signs);
        volelem = abs(det([loccoords(:,2)-loccoords(:,1),loccoords(:,3)-loccoords(:,1)]))/2;
        diffcoords = loccoords(:)*ones(1,3) - repmat(loccoords,3,1);
        patlocnodes = reshape(diffcoords*voledgeswsigns*guref(globedgeinds),2,3)'/(2*volelem);
    
        % write local nodal values into global vector
        guref_x(elem,:) = patlocnodes(:,1)';
        guref_y(elem,:) = patlocnodes(:,2)';
    end % end
    guref_x = reshape(guref_x',[],1);
    guref_y = reshape(guref_y',[],1);
    
    % dg0 and dg1 mass matrices
    Mdg0 = Mref;
    Mdg1 = kron(spdiags(diag(Mref),0,Tref.nelems,Tref.nelems),[2 1 1;1 2 1;1 1 2]./12);
        
    %% compute errors against the reference solution
    for level = 0:maxlevel
        fprintf('computing for h = 2^{-%d}\n',level);
        % mesh
        T = refineMesh(getMesh('scaledSquareSymFKrot',L),level);
        [~,P0,P1] = refineMesh(T,nref-level);
        % assemble RT0 matrices
        [edges,nodes2edges] = getEdgeProperties(T);
        [B,C,M] = assembleRT0(T);
    
        % evaluate coarse potential
        potential_eval = potential(computeMids(T));
        
        % nonlinear eigenvalue solver
        % initial iteration
        mids = computeMids(T);
        u0 = (1-mids(:,1).^2/L^2).*(1-mids(:,2).^2/L^2);
        u0 = u0/sqrt(u0'*(M*u0));
        
        % [u,lambda,energy,normres] = Amethod(u0,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs);
        [u,lambda,energy,normres] = Jmethod(u0,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs);
    
        % compute gradient of u
        [gu,~,~,~] = pcg(B,-C'*u);
    
        % reconstruct flux from vector gu
        gu_x = zeros(T.nelems,3);
        gu_y = zeros(T.nelems,3);
        for elem = 1:T.nelems
            % evaluate RT0 basis functions at local nodes
            globnodes = T.elems(elem,:);
            globedgeinds = full(nodes2edges(sub2ind([T.nnodes,T.nnodes],globnodes([2 1 1]),globnodes([3 3 2]))));
            loccoords = T.coords(globnodes,:)';
            voledges = [norm(loccoords(:,2)-loccoords(:,3),2), ...
                        norm(loccoords(:,1)-loccoords(:,3),2),...
                        norm(loccoords(:,1)-loccoords(:,2),2)];
            signs = -1 + 2*(edges(globedgeinds,3) == elem)';
            voledgeswsigns = diag(voledges.*signs);
            volelem = abs(det([loccoords(:,2)-loccoords(:,1),loccoords(:,3)-loccoords(:,1)]))/2;
            diffcoords = loccoords(:)*ones(1,3) - repmat(loccoords,3,1);
            patlocnodes = reshape(diffcoords*voledgeswsigns*gu(globedgeinds),2,3)'/(2*volelem);
        
            % write local nodal values into global vector
            gu_x(elem,:) = patlocnodes(:,1)';
            gu_y(elem,:) = patlocnodes(:,2)';
        end % end
        gu_x = reshape(gu_x',[],1);
        gu_y = reshape(gu_y',[],1);
    
        % compute errors
        % u
        u_diff = uref-P0*u;
        u_err(indkappa,level+1) = sqrt(u_diff'*(Mdg0*u_diff))/sqrt(uref'*(Mdg0*uref));
        % gu
        gu_x_diff = guref_x - P1*gu_x;
        gu_y_diff = guref_y - P1*gu_y;
        gu_err(indkappa,level+1) = sqrt(gu_x_diff'*(Mdg1*gu_x_diff) + gu_y_diff'*(Mdg1*gu_y_diff))/sqrt(guref_x'*(Mdg1*guref_x) + guref_y'*(Mdg1*guref_y));
        % lambda
        lambda_err(indkappa,level+1) = abs(lambdaref(end) - lambda(end))/abs(lambdaref(end));
        % energy
        energy_err(indkappa,level+1) = abs(energyref(end) - energy(end))/abs(energyref(end));
    end % for
    % plot solution
    f0 = figure('position',[100,100,500,500]);
    plotDiscFun(Tref,uref,'edgecolor','none');
    axis('square');
    xticks(linspace(-L,L,5)); xlim([-L,L]);
    yticks(linspace(-L,L,5)); ylim([-L,L]);
    colorbar
        
    exportgraphics(f0,['plots/groundstate_harmpot_kappa' num2str(kappa) '.png']);
end % for

%% plot the results 
% convergence plot for u and gu
f1 = figure('position',[100,100,500,500]);
hs = 2.^-(0:maxlevel);
loglog(hs,4*hs,':k');
hold on;
legendcell = {'slope 1'};
for indkappa = 1:length(kappas)
    kappa = kappas(indkappa);
    loglog(hs,u_err(indkappa,:),'color',[0 0.4470 0.7410],'marker','x','markersize',indkappa+4);
    legendcell = [legendcell,['$\|u^\mathrm{ref}-u_h\|_{L^2}/\|u^\mathrm{ref}\|_{L^2},\;\kappa = ' num2str(kappa) '$']];
end % for
for indkappa = 1:length(kappas)
    kappa = kappas(indkappa);
    loglog(hs,gu_err(indkappa,:),'color',[0.8500 0.3250 0.0980],'marker','o','markersize',indkappa+4);
    legendcell = [legendcell,['$\|\nabla u-\sigma_h\|_{L^2}/\|\nabla u^\mathrm{ref}\|_{L^2},\;\kappa = ' num2str(kappa) '$']];
end % for
legend(legendcell,'interpreter','latex');
xlabel('$h$','interpreter','latex');
xticks(flip(hs));
xticklabels({'2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'})

exportgraphics(f1,'plots/convergencekappa_ugu.png');

% convergence plot for lambda and energy
f2 = figure('position',[100,100,500,500]);
loglog(hs,.5*hs.^2,':k');
hold on;
legendcell = {'slope 2'};
for indkappa = 1:length(kappas)
    kappa = kappas(indkappa);
    loglog(hs,energy_err(indkappa,:),'color',[0 0.4470 0.7410],'marker','x','markersize',indkappa+4);
    legendcell = [legendcell,['$|E^\mathrm{ref}-E_h|/E^\mathrm{ref},\;\kappa = ' num2str(kappa) '$']];
end % for
for indkappa = 1:length(kappas)
    kappa = kappas(indkappa);
    loglog(hs,lambda_err(indkappa,:),'color',[0.8500 0.3250 0.0980],'marker','o','markersize',indkappa+4);
    legendcell = [legendcell,['$|\lambda-\lambda_h|/\lambda^\mathrm{ref},\;\kappa = ' num2str(kappa) '$']];
end % for
legend(legendcell,'interpreter','latex');
xlabel('$h$','interpreter','latex');
xticks(flip(hs));
xticklabels({'2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'})

exportgraphics(f2,'plots/convergencekappa_energylambda.png');
