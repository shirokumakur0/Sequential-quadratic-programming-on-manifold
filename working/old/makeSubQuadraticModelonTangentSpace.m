function subproblem = makeSubQuadraticModelonTangentSpace(problem, xCur, mus, lambdas)

    N = tangentspacefactory(problem.M, xCur);
    subproblem.M = N;
    % costLag = @(X) costLagrangian(X, problem, mus, lambdas); % though we don't use here.
    gradLag = @(X) gradLagrangian(X, problem, mus, lambdas); % in the tangent space
    hessLag = @(X, d) hessLagrangian(X, d, problem, mus, lambdas); % in the tangent space
    
    subproblem.cost = @(d) .5 * subproblem.M.inner(xCur, hessLag(xCur, d), d)...
                           + subproblem.M.inner(xCur, gradLag(xCur), d);
    subproblem.grad = @(d) hessLag(xCur, d) + gradLag(xCur);
    % Note that grad above is already 'Riemannian' version since gradLag and
    % hessLag return the tangent vector, respectively; see their programs for
    % more details.
    
    if problem.condet.has_ineq_cost
        n_ineq_constraint = problem.condet.n_ineq_constraint_cost;
        subproblem.ineq_constraint_cost = cell(n_ineq_constraint, 1);
        subproblem.ineq_constraint_grad = cell(n_ineq_constraint, 1);
        for numineq = 1: n_ineq_constraint
            costhandle = problem.ineq_constraint_cost{numineq};
            gradhandle = problem.ineq_constraint_grad{numineq};
            subproblem.ineq_constraint_cost{numineq} = @(d) subproblem.M.inner(xCur,...
            gradhandle(xCur), d) + costhandle(xCur);
            subproblem.ineq_constraint_grad{numineq} = @(d) gradhandle(xCur);
            % the gradients above are all egrad vectors! watch out!
        end
    end
    
    if problem.condet.has_eq_cost
        n_eq_constraint = problem.condet.n_eq_constraint_cost;
        subproblem.eq_constraint_cost = cell(n_eq_constraint, 1);
        subproblem.eq_constraint_grad = cell(n_eq_constraint, 1);
        for numeq = 1: n_eq_constraint
            costhandle = problem.eq_constraint_cost{numeq};
            gradhandle = problem.eq_constraint_grad{numeq};
            subproblem.eq_constraint_innercost{numeq} = @(d) subproblem.M.inner(xCur, gradhandle(xCur), d);
            subproblem.eq_constraint_constcost{numeq} = costhandle(xCur);
            subproblem.eq_constraint_cost{numeq} = @(d) subproblem.M.inner(xCur,...
                gradhandle(xCur), d) + costhandle(xCur);
            subproblem.eq_constraint_grad{numeq} = @(d) gradhandle(xCur);
        end
    end
    % the gradients above are all egrad vectors! watch out!

    subproblem.condet = constraintsdetail(subproblem);

end