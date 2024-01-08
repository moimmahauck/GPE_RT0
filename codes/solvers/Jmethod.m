function [u,lambda,energy,normres] = Jmethod(u,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs)
    lambda = zeros(1,maxiter);
    energy = zeros(1,maxiter);
    normres = zeros(1,maxiter);
    
    alpha = zeros(1,3);
    beta = zeros(1,5);
    zeta = zeros(1,3);
    
    V_M = spdiags(M,0)';
    M_inv = spdiags(V_M(:).^(-1),0,size(M,1),size(M,1));
    M_nl_inv = @(u,sigma) spdiags((V_M'.*(potential_eval + 3*kappa*u.^2 + sigma)).^(-1),0,size(M,1),size(M,1));
    x = @(u) -2*kappa*V_M'.*u.^3;
    y = @(u) V_M'.*u;
    [BinvCtu,~,~,~] = pcg(B,C'*u);
    
    for iter = 1:maxiter
        % turn on shifting if residual is smaller than 1e-2
        if iter == 1 || normres(iter-1) > 1e-2
            sigma = 0;
        else
            sigma = -lambda(iter-1);
        end % if
        M_nl_inv_eval = M_nl_inv(u,sigma);
        % woodbury matrix identity used to simplify the iteration
        S = B + C'*(M_nl_inv_eval*C);
        tmp = M_nl_inv_eval*[x(u) M*u] - M_nl_inv_eval*C*(S\(C'*(M_nl_inv_eval*[x(u) M*u])));
        tmpx = tmp(:,1); tmpMu = tmp(:,2);
        JsinvMu = tmpMu - tmpx*(y(u)'*tmpMu)/(1+y(u)'*tmpx);
        gamma = (JsinvMu'*(M*u)).^(-1);
        if iter == 1 || normres(iter-1) > 1e-1
            % perform golden section search 
            [BinvCtJsinvMu,~,~,~] = pcg(B,C'*JsinvMu);
            % compute values of alpha, beta, and zeta
            alpha(1) = u'*(C*BinvCtu) + V_M*(potential_eval.*u.^2);
            alpha(2) = 2*gamma*(u'*(C*BinvCtJsinvMu) + V_M*(potential_eval.*u.*JsinvMu));
            alpha(3) = gamma.^2*(JsinvMu'*(C*BinvCtJsinvMu) + V_M*(potential_eval.*JsinvMu.^2));
            beta(1) = .5*kappa*V_M*(u.^4);
            beta(2) = 2*kappa*gamma*V_M*(u.^3.*JsinvMu);
            beta(3) = 3*kappa*gamma^2*V_M*(u.^2.*JsinvMu.^2);
            beta(4) = 2*kappa*gamma^3*V_M*(u.*JsinvMu.^3);
            beta(5) = .5*kappa*gamma^4*V_M*(JsinvMu.^4);
            zeta(1) = u'*(M*u);
            zeta(2) = 2*gamma*u'*(M*JsinvMu);
            zeta(3) = gamma^2*JsinvMu'*(M*JsinvMu);
            % energy in dependence of tau
            f = @(tau) computeEnergy(tau,alpha,beta,zeta);
            % golden section search in the interval [0,2]
            tau = goldenSearch(f,0,2,tol_gs,maxiter_gs);
        else 
            tau = 1;
        end % if
        % new iterate
        u = (1-tau)*u + tau*gamma*JsinvMu;
        u = u/sqrt(u'*(M*u));
    
        % compute residual, eigenvalue and energy
        [BinvCtu,~,~,~] = pcg(B,C'*u);
        energy(iter) = .5*u'*(C*BinvCtu) + .5*V_M*(potential_eval.*u.^2) + .25*kappa*V_M*(u.^4);
        lambda(iter) = 2*energy(iter) + .5*kappa*V_M*(u.^4);
        % residual in L2'-norm
        res = C*BinvCtu + V_M'.*(potential_eval + kappa*u.^2).*u - lambda(iter)*M*u;
        normres(iter) = res'*(M_inv*res);
        fprintf('iteration: %d, energy: %2.8f, residual: %d\n',iter,energy(iter),normres(iter));
        if iter > 1 && abs(energy(iter)-energy(iter-1))/energy(iter) <= tol
            break;
        end % if
    end % for
    lambda = lambda(1:iter); energy = energy(1:iter); normres = normres(1:iter);
end % for