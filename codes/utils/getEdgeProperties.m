function [edges,nodes2edges] = getEdgeProperties(T)
    % compute edges and associate elements to edges
    edgestmp = [T.elems(:,[1,2]); T.elems(:,[1,3]); T.elems(:,[2,3])];
    edgestmp = [sort(edgestmp,2) repmat(1:T.nelems,1,3)'];
    edgestmp = sortrows(edgestmp,[1 2]);
    % eliminte doubplicate edges
    [edgesunqiue,ind1,ind2] = unique(edgestmp(:,1:2),'rows','stable');
    [~,tmp] = unique(ind2,'stable');
    inddouplicates = setdiff(1:numel(ind2),tmp);
    % add the adjacent element of doubplicate edges
    edges = [edgesunqiue edgestmp(ind1,3) zeros(size(edgesunqiue,1),1)];
    edges(ind2(inddouplicates),4) = edgestmp(inddouplicates,3);

    % compute matrix that associates pairs of nodes to edges 
    nodes2edges = sparse(edges(:,1),edges(:,2),1:size(edges,1),T.nnodes,T.nnodes);
    nodes2edges = nodes2edges + nodes2edges';
end % function