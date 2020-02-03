function H = innerHess2matrix(quadcost, basis)
    n = numel(basis);
    H = zeros(n);
    for i=1:n
        H(i,i) = quadcost(basis{i}, basis{i});
        for j = (i+1):n
            H(i,j) = quadcost(basis{i}, basis{j});
            H(j,i) = H(i,j);
        end
    end
end