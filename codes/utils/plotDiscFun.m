function plotDiscFun(T,x,varargin)
    % P0 input
    if size(x,2) == 1
        x = repmat(x,1,3);
    end % if

    % coordinates for all nodes of sides
    X1 = T.coords(T.elems(:,1), 1);
    X2 = T.coords(T.elems(:,2), 1);
    X3 = T.coords(T.elems(:,3), 1);
    Y1 = T.coords(T.elems(:,1), 2);
    Y2 = T.coords(T.elems(:,2), 2);
    Y3 = T.coords(T.elems(:,3), 2);
    
    % quads for the function in patch-style
    valX = [X1';X2';X3'];
    valY = [Y1';Y2';Y3'];
    
    % nodal values at each element
    valZ = [x(:,1)';x(:,2)';x(:,3)'];

    patch(valX,valY,valZ,valZ,varargin{:});
end