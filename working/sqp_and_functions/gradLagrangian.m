function gradLag = gradLagrangian(x, problem0, mus, lambdas)
    gradLag = getGradient(problem0, x);
    if problem0.condet.has_ineq_cost
        for numineq = 1: problem0.condet.n_ineq_constraint_cost
            gradhandle = problem0.ineq_constraint_grad{numineq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
            gradLag = problem0.M.lincomb(x, 1, gradLag, mus(numineq), constraint_grad);
        end
    end

    if problem0.condet.has_eq_cost
        for numeq = 1:problem0.condet.n_eq_constraint_cost
            gradhandle = problem0.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
            gradLag = problem0.M.lincomb(x, 1, gradLag, lambdas(numeq), constraint_grad);
        end
    end
end