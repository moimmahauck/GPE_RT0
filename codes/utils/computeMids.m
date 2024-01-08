function mids = computeMids(T)
    mids = (T.coords(T.elems(:,1),:) + T.coords(T.elems(:,2),:) + T.coords(T.elems(:,3),:))/3;
end % funciton