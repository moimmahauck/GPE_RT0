function [u,lambda,energy,normres] = Amethod(u,B,C,M,potential_eval,kappa,tol,maxiter,tol_gs,maxiter_gs)
    lambda = zeros(1,maxiter);
    energy = zeros(1,maxiter);
    normres = zeros(1,maxiter);
    
    alpha = zeros(1,3);
    beta = zeros(1,5);
    zeta = zeros(1,3);

    V_M = spdiags(M,0)';
    M_inv = spdiags(V_M(:).^(-1),0,size(M,1),size(M,1));
    M_nl = @(u) spdiags(V_M'.*(potential_eval + kappa*u.^2),0,size(M,1),size(M,1));
    M_nl_inv = @(u) spdiags((V_M'.*(potential_eval + kappa*u.^2)).^(-1),0,size(M,1),size(M,1));
    
    [BinvCtu,~,~,~] = pcg(B,C'*u);
    
    for iter = 1:maxiter  
        M_nl_inv_eval = M_nl_inv(u);
        % woodbury matrix identity used to simplify the iteration
        S = B + C'*(M_nl_inv_eval*C);
        Gu = M_nl_inv_eval*M*u - M_nl_inv_eval*C*(S\(C'*(M_nl_inv_eval*M*u)));
        [BinvCtGu,~,~,~] = pcg(B,C'*Gu);
        gamma = (u'*(M*u))/(Gu'*(C*BinvCtGu) + Gu'*(M_nl(u)*Gu));
        % compute values of alpha, beta, and zeta
        alpha(1) = u'*(C*BinvCtu) + V_M*(potential_eval.*u.^2);
        alpha(2) = 2*gamma*(u'*(C*BinvCtGu) + V_M*(potential_eval.*u.*Gu));
        alpha(3) = gamma.^2*(Gu'*(C*BinvCtGu) + V_M*(potential_eval.*Gu.^2));
        beta(1) = .5*kappa*V_M*(u.^4);
        beta(2) = 2*kappa*gamma*V_M*(u.^3.*Gu);
        beta(3) = 3*kappa*gamma^2*V_M*(u.^2.*Gu.^2);
        beta(4) = 2*kappa*gamma^3*V_M*(u.*Gu.^3);
        beta(5) = .5*kappa*gamma^4*V_M*(Gu.^4);
        zeta(1) = u'*(M*u);
        zeta(2) = 2*gamma*u'*(M*Gu);
        zeta(3) = gamma^2*Gu'*(M*Gu);
        % energy in dependence of tau
        f = @(tau) computeEnergy(tau,alpha,beta,zeta);
        % golden section search in the interval [0,2]
        tau = goldenSearch(f,0,2,tol_gs,maxiter_gs);
        % new iterate
        u = (1-tau)*u + tau*gamma*Gu;
        u = u/sqrt(u'*(M*u));
        % compute lambda and energy
        energy(iter) = f(tau);
        lambda(iter) = 2*energy(iter) + .5*kappa*V_M*(u.^4);
        % residual in L2'-norm
        [BinvCtu,~,~,~] = pcg(B,C'*u);
        res = C*BinvCtu + M_nl(u)*u - lambda(iter)*M*u;
        normres(iter) = res'*(M_inv*res);
        fprintf('iteration: %d, energy: %2.8f, residual: %d\n',iter,energy(iter),normres(iter));
        if iter > 1 && abs(energy(iter)-energy(iter-1))/energy(iter) <= tol
            break;
        end % if
    end % for
    lambda = lambda(1:iter); energy = energy(1:iter); normres = normres(1:iter);
end % function