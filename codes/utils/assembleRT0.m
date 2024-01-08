function [B,C,M] = assembleRT0(T)
    [edges,nodes2edges] = getEdgeProperties(T);
    nedges = size(edges,1);

    % matrices for storing local indices and local matrices
    I_B = zeros(3,3*T.nelems);
    J_B = zeros(3,3*T.nelems);
    V_B = zeros(3,3*T.nelems);
    
    I_C = zeros(3,T.nelems);
    J_C = zeros(3,T.nelems);
    V_C = zeros(3,T.nelems);
    
    V_M = zeros(1,T.nelems);
    
    % assemble local matrices
    factors = [2 0 1 0 1 0;0 2 0 1 0 1;1 0 2 0 1 0;0 1 0 2 0 1;1 0 1 0 2 0;0 1 0 1 0 2];
    for elem = 1:T.nelems
        globnodes = T.elems(elem,:);
        globedgeinds = full(nodes2edges(sub2ind([T.nnodes,T.nnodes],globnodes([2 1 1]),globnodes([3 3 2]))));
        loccoords = T.coords(globnodes,:)';
        diffcoords = loccoords(:)*ones(1,3) - repmat(loccoords,3,1);
        voledges = [norm(loccoords(:,2)-loccoords(:,3),2), ...
                    norm(loccoords(:,1)-loccoords(:,3),2),...
                    norm(loccoords(:,1)-loccoords(:,2),2)];
        signs = -1 + 2*(edges(globedgeinds,3) == elem)';
        voledgeswsigns = diag(voledges.*signs);
        volelem = abs(det([loccoords(:,2)-loccoords(:,1),loccoords(:,3)-loccoords(:,1)]))/2;
    
        V_B(:,3*(elem-1)+1:3*elem) = voledgeswsigns*diffcoords'*factors*diffcoords*voledgeswsigns/(48*volelem);
        I_B(:,3*(elem-1)+1:3*elem) = repmat(globedgeinds',1,3);
        J_B(:,3*(elem-1)+1:3*elem) = repmat(globedgeinds,3,1);
    
        V_C(:,elem) = (voledges.*signs)';
        I_C(:,elem) = elem;
        J_C(:,elem) = globedgeinds';
    
        V_M(elem) = volelem;
    end % for
    
    % assemble matrices B and C    
    B = sparse(I_B(:),J_B(:),V_B(:),nedges,nedges); 
    C = sparse(I_C(:),J_C(:),V_C(:),T.nelems,nedges);
    M = spdiags(V_M(:),0,T.nelems,T.nelems);
end % for
