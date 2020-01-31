function canonicalbasis = makeCanonicalBasis(problem0)
    % Create basis of tangent spaces whose dimension is eqaul to
    % problem0.M, which is the original manifold. The shape of cells and
    % each component, that is a basis, follow the original manifold. In
    % other words, if points on a manifold is a vector, then a cell and basis
    % are also vectors. If points on a manifold is matrix, then so are a
    % cell and basis.
    
    % Here we assume that the tangent spaces
    % of M have the same basis as that of M. At least, this seems to be true
    % when M is embedded in some Euclidean space unless factories in Manopt
    % haven't changed after Jan. 27, 2020, i.e., both M and the tangent spaces
    % use the canonical basis. Yet, we should consider the consistency
    % between M and the tangent sapces when applying this method to some M,
    % respectively. 
    % This part is an uncomplete if we launch this algorithm, officialy.
    
    canonicalbasis = cell(size(problem0.M.zerovec()))
    n = numel(problem0.M.zerovec())
    for k = 1 : n
        vec = problem0.M.zerovec();
        vec(k) = 1;
        canonicalbasis{k} = vec;
    end
end