function [T,P0,P1] = refine(T)
    edges = [T.elems(:,[1,2]); T.elems(:,[1,3]); T.elems(:,[2,3])];
    edges = sort(edges,2);
    d2p = sparse(edges(:,1),edges(:,2),1,T.nnodes,T.nnodes);
    [edges1,edges2] = find(d2p);
    edges = [edges1,edges2];
    nedges = size(edges,1);

    % new coordinates
    newcoords = 0.5.*(T.coords(edges(:,1),:)+T.coords(edges(:,2),:)); 
    T.coords = [T.coords; newcoords]; 

    % new elements
    edgemarker = (T.nnodes+1:T.nnodes+nedges)';
    d2p = sparse(edges,edges(:,[2,1]),(1:nedges)'*[1 1],T.nnodes,T.nnodes);
    elemstmp = zeros(T.nelems,6);
    elemstmp(:,1:3) = T.elems; 
    elemstmp(:,4) = edgemarker(d2p(T.elems(:,1)+T.nnodes*(T.elems(:,2)-1)));
    elemstmp(:,5) = edgemarker(d2p(T.elems(:,2)+T.nnodes*(T.elems(:,3)-1)));
    elemstmp(:,6) = edgemarker(d2p(T.elems(:,3)+T.nnodes*(T.elems(:,1)-1)));
    T.elems = [elemstmp(:,[1,4,6])
               elemstmp(:,[4,2,5])
               elemstmp(:,[6,5,3])
               elemstmp(:,[4,5,6])];

    % P0 prolongation matrix
    new2old = repmat(1:T.nelems,1,4)';
    P0 = sparse((1:length(new2old))',new2old,1,length(new2old),T.nelems);
    
    % P1 discontinuous porolongation matrix
    P1 = [kron(speye(T.nelems),[1 0 0;.5 .5 0;.5 0 .5])
          kron(speye(T.nelems),[.5 .5 0;0 1 0;0 .5 .5])
          kron(speye(T.nelems),[.5 0 .5;0 .5 .5;0 0 1])
          kron(speye(T.nelems),[.5 .5 0;0 .5 .5;.5 0 .5])];
    
    % update number of nodes and edges
    T.nnodes = size(T.coords,1);
    T.nelems = size(T.elems,1);
end % function