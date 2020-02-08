function [ineqconst_gradmat, ineqconst_costvec, ...
     eqconst_gradmat, eqconst_costvec] = gradConstraintMatrix(xCur, problem0, basis)
    
    ineqconst_gradmat = [];
    ineqconst_costvec = [];
    eqconst_gradmat = [];
    eqconst_costvec = [];
    if problem0.condet.has_ineq_cost
        row = numel(basis);
        col = problem0.condet.n_ineq_constraint_cost;
        for numineq = 1: problem0.condet.n_ineq_constraint_cost
            gradhandle = problem0.ineq_constraint_grad{numineq}; % can be refactored
            constraint_grad = gradhandle(xCur);
            constraint_grad = problem0.M.egrad2rgrad(xCur, constraint_grad);
            constraint_gradvec = gradMetricVectorize(xCur, constraint_grad, problem0, basis);
            constraint_gradvec = transpose(constraint_gradvec);
            ineqconst_gradmat = vertcat(ineqconst_gradmat, constraint_gradvec);            
            costhandle = problem0.ineq_constraint_cost{numineq};  % can be refactored
            ineqconst_costvec(numineq) = costhandle(xCur);           
        end
    end

    if problem0.condet.has_eq_cost
        row = numel(basis);
        col = problem0.condet.n_eq_constraint_cost;
        for numeq = 1 : problem0.condet.n_eq_constraint_cost
            gradhandle = problem0.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(xCur);
            constraint_grad = problem0.M.egrad2rgrad(xCur, constraint_grad);
            constraint_gradvec = gradMetricVectorize(xCur, constraint_grad, problem0, basis);
            constraint_gradvec = transpose(constraint_gradvec);
            eqconst_gradmat = vertcat(eqconst_gradmat, constraint_gradvec);
            costhandle = problem0.eq_constraint_cost{numeq};
            eqconst_costvec(numeq) = costhandle(xCur);
        end
    end
end