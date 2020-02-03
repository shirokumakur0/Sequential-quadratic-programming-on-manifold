function egradLag = egradLagrangian(x, problem0, mus, lambdas)
    egradLag = problem0.egrad(x);
    if problem0.condet.has_ineq_cost
        for numineq = 1: problem0.condet.n_ineq_constraint_cost
            gradhandle = problem0.ineq_constraint_grad{numineq};
            constraint_grad = gradhandle(x);
            egradLag = problem0.M.lincomb(x, 1, egradLag, mus(numineq), constraint_grad);
        end
    end

    if problem0.condet.has_eq_cost
        for numeq = 1:problem0.condet.n_eq_constraint_cost
            gradhandle = problem0.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(x);
            egradLag = problem0.M.lincomb(x, 1, egradLag, lambdas(numeq), constraint_grad);
        end
    end
end