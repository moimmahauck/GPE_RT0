function T = getMesh(geom,factor)
    switch geom
        case 'unitSquareFK'
            elems = importdata('elems0.dat');
            coords = importdata('coords0.dat');
        case 'scaledSquareFK'
            elems = importdata('elems0.dat');
            coords = importdata('coords0.dat');
            coords = factor*(2*coords-ones(1,2));
        case 'scaledSquareSym'
            elems = importdata('elems1.dat');
            coords = importdata('coords1.dat');
            coords = factor*coords;
        case 'scaledSquareSymFKrot'
            elems = importdata('elems2.dat');
            coords = importdata('coords2.dat');
            coords = factor*coords;
        otherwise
            error('geometry not found.');
    end % switch

    % create mesh
    T = struct('coords',coords,'elems',elems,'nnodes',size(coords,1),'nelems',size(elems,1));
end % function