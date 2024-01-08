function [T,P0,P1] = refineMesh(T,nref)
    P0 = speye(T.nelems,T.nelems);
    P1 = speye(3*T.nelems,3*T.nelems);
    for k=1:nref
        [T,P0lvl,P1lvl] = refine(T);
        P0 = P0lvl*P0;
        P1 = P1lvl*P1;
    end % for
end % function
