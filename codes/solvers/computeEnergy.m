function energy = computeEnergy(tau,alpha,beta,zeta)
    s = 0;
    for i = 0:2
        j = 2 - i;
        s = s + (1-tau)^i*tau^j*zeta(j+1);
    end % for
    s = 1/sqrt(s);

    energy1 = 0;
    for i = 0:2
        j = 2 - i;
        energy1 = energy1 + s^2*(1-tau)^i*tau^j*alpha(j+1);
    end % for

    energy2 = 0;
    for i = 0:4
        j = 4 - i;
        energy2 = energy2 + s^4*(1-tau)^i*tau^j*beta(j+1);
    end % for
    energy = .5*energy1 + .5*energy2;
end % funciton