function hessLag = hessLagrangian(x, dir, problem0, mus, lambdas)
    hessLag = getHessian(problem0, x, dir);
    if problem0.condet.has_ineq_cost
        for numineq = 1 : problem0.condet.n_ineq_constraint_cost
            gradhandle = problem0.ineq_constraint_grad{numineq};
            constraint_egrad = gradhandle(x); % to be refactored
            hesshandle = problem0.ineq_constraint_hess{numineq};
            constraint_ehess = hesshandle(x, dir);
            constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                     constraint_ehess, dir);
            hessLag = problem0.M.lincomb(x, 1, hessLag,...
                mus(numineq), constraint_hess);
        end
    end
    if problem0.condet.has_eq_cost
        for numeq = 1 : problem0.condet.n_eq_constraint_cost
            gradhandle = problem0.eq_constraint_grad{numeq};
            constraint_egrad = gradhandle(x);
            hesshandle = problem0.eq_constraint_hess{numeq};
            constraint_ehess = hesshandle(x, dir);
            constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                    constraint_ehess, dir);
            hessLag = problem0.M.lincomb(x, 1, hessLag,...
                lambdas(numeq), constraint_hess);
        end
    end
end