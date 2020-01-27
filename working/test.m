% n = M.dim()
% basis = cell(n,1)
% for k=1:n
%    vec = zeros(n,1);
%    vec(k) = 1;
%    basis{k} = vec
% end


function [ineqconst_gradmat, eqconst_gradmat] = gradConstraintMatrix(xCur, problem0, basis)
    if condet.has_ineq_cost
        row = numel(basis);
        col = condet.n_ineq_constraint_cost;
        
        for numineq = 1: condet.n_ineq_constraint_cost
            % costhandle = problem.ineq_constraint_cost{numineq};
            % cost_numineq = costhandle(x);
            gradhandle = problem0.ineq_constraint_grad{numineq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
            val = problem0.M.lincomb(x, 1, val, mus(numineq), constraint_grad);
        end
    end

    if condet.has_eq_cost
        for numeq = 1:condet.n_eq_constraint_cost
            % costhandle = problem.eq_constraint_cost{numeq};
            % cost_numeq = costhandle(x);
            gradhandle = problem0.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
            val = problem0.M.lincomb(x, 1, val, lambdas(numeq), constraint_grad);
        end
    end
end