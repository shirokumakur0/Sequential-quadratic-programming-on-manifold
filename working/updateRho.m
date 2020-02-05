function rho = updateRho(rho, Lagmultipliers, problem0)

    newacc = 0;
    for iterineq = 1 : problem0.condet.n_ineq_constraint_cost
        newacc = max(newacc, Lagmultipliers.ineqlin(iterineq));
    end
    for itereq = 1 : problem0.condet.n_eq_constraint_cost
        newacc = max(newacc, abs(Lagmultipliers.eqlin(itereq)));
    end

    if rho < newacc
       rho = newacc;
    end
    
end