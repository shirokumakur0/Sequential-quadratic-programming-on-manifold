function gradmetvec = gradMetricVectorize(x, grad, problem, basis)
    n = numel(basis); 
    gradmetvec = zeros(n,1);
    for i = 1 : n
            gradmetvec(i) = problem.M.inner(x, grad, basis{i});
    end
end