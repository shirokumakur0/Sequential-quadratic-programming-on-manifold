function qpinfo = makeQPInfo(problem, xCur, mus, lambdas, options)
    costLag = @(X) costLagrangian(X, problem, mus, lambdas); % though we don't use here.
    gradLag = @(X) gradLagrangian(X, problem, mus, lambdas); % in the tangent space
    hessLag = @(X, d) hessLagrangian(X, d, problem, mus, lambdas); % in the tangent space

    auxproblem.M = problem.M;
    auxproblem.cost = costLag;
    auxproblem.grad = gradLag;
    auxproblem.hess = hessLag;
    
    %If we'll regularize Hessian in the way of 'mineigval_manopt',
    %calculate the minimum eigenvalue here.
    if strcmp(options.trimhessian, "mineigval_manopt")
        [~ ,qpinfo.mineigval_manopt] = hessianextreme(auxproblem, xCur);
    end
    
    % make H and basis
    [H,basis] = hessianmatrix(auxproblem, xCur);
    qpinfo.H = H;
    qpinfo.basis = basis;
    qpinfo.n = numel(basis);
    
    % make f
    f = zeros(qpinfo.n, 1);
    for i=1:qpinfo.n
        f(i) = auxproblem.M.inner(xCur, auxproblem.grad(xCur), basis{i});
    end
    qpinfo.f = f;
    
    % make inequality constraints
    if problem.condet.has_ineq_cost
        row = problem.condet.n_ineq_constraint_cost;
        col = numel(basis);
        A = zeros(row, col);
        b = zeros(row, 1);
        for i = 1:row
            costhandle = problem.ineq_constraint_cost{i};
            b(i) = - costhandle(xCur);
            gradhandle = problem.ineq_constraint_grad{i};
            constraint_egrad = gradhandle(xCur);
            constraint_grad = problem.M.egrad2rgrad(xCur, constraint_egrad);
            for j = 1:col
                base = basis{j};
                A(i,j) = problem.M.inner(xCur, constraint_grad, base);
            end
        end
    else
        A = [];
        b = [];
    end
    qpinfo.A = A;
    qpinfo.b = b;
    
    % make equality constraints
    if problem.condet.has_eq_cost
        row = problem.condet.n_eq_constraint_cost;
        col = numel(basis);
        Aeq = zeros(row, col);
        beq = zeros(row, 1);
        for i = 1:row
            costhandle = problem.eq_constraint_cost{i};
            beq(i) = - costhandle(xCur);
            gradhandle = problem.eq_constraint_grad{i};
            constraint_egrad = gradhandle(xCur);
            constraint_grad = problem.M.egrad2rgrad(xCur, constraint_egrad);
            for j = 1:col
                base = basis{j};
                Aeq(i,j) = problem.M.inner(xCur, constraint_grad, base);
            end
        end
    else
        Aeq = [];
        beq = [];
    end
    qpinfo.Aeq = Aeq;
    qpinfo.beq = beq;
end