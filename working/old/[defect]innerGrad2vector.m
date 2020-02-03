function f = innerGrad2vector(lincost, basis)
    n = numel(basis);
    f = zeros(n, 1);
    for i=1:n
        f(i) = lincost(basis{i});
    end
end