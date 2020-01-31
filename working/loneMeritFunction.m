function val = loneMeritFunction(problem0, x, rho)

    val = getCost(problem0, x);
    if problem0.condet.has_ineq_cost
        for numineq = 1: problem0.condet.n_ineq_constraint_cost
            costhandle = problem0.ineq_constraint_cost{numineq};
            cost_numineq = costhandle(x);
            val = val + rho * max(0, cost_numineq);
        end
    end

    if problem0.condet.has_eq_cost
        for numeq = 1: problem0.condet.n_eq_constraint_cost
            costhandle = problem0.eq_constraint_cost{numeq};
            cost_numeq = costhandle(x);
            val = val + rho * abs(cost_numeq);
        end
    end
end
