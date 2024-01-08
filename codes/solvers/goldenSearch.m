function x = goldenSearch(f,a,b,tol,maxiter)
    gr = .5*(sqrt(5) - 1);
    iter = 1;
    while true
        x1 = a + gr*(b - a); x2 = b - gr*(b - a);
        if f(x1) >= f(x2)
            b = x1;
        else
            a = x2;
        end % if
        if iter > maxiter || abs(x1 -x2) < tol
            break;
        end % if
        iter = iter + 1;
    end % while
    x = .5*x1 + .5*x2;
end % function