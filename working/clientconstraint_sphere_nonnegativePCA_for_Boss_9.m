function data = clientconstraint_sphere_nonnegativePCA(X, rankY, methodoptions, specifier)
%Y is a square matrix
%rank is smaller than columns of Y
[N, ~] = size(X);

data = NaN(3, 5);

M = spherefactory(N, rankY);
problem.M = M;
problem.cost = @(u) costfun(u);
problem.egrad = @(u) gradfun(u);
problem.ehess = @(u,d) hessfun(u,d);
x0 = M.rand();

%     DEBUG only
%     checkgradient(problem);

%-------------------------Set-up Constraints-----------------------
constraints_cost = cell(1, N*rankY);
for row = 1: N
    for col = 1: rankY
        constraints_cost{(col-1)*N + row} = @(Y) -Y(row, col);
    end
end

constraints_grad = cell(1, N * rankY);
for row = 1: N
    for col = 1: rankY
        constraintgrad = zeros(N, rankY);
        constraintgrad(row, col) = -1;
        constraints_grad{(col-1)*N + row} = @(U) constraintgrad;
    end
end

constraints_hess = cell(1, N * rankY);
for row = 1: N
    for col = 1: rankY
        constraints_hess{(col-1)*N + row} = @(U, D) 0;
    end
end

problem.ineq_constraint_cost = constraints_cost;
problem.ineq_constraint_grad = constraints_grad;
problem.ineq_constraint_hess = constraints_hess;
%     Debug Only
%     checkconstraints(problem)

condet = constraintsdetail(problem);

%     ------------------------- Solving ---------------------------
    options = methodoptions;
        
    fprintf('Starting SQP\n');
    
    [x, xcost, info] = sqp(problem, x0, options); % #ok<ASGLU>   

     
%------------------------sub functions-----------     
     
     function stop = outfun(x, optimValues, state)
        stop = false;
        if toc(timetic) > methodoptions.maxtime
            stop = true;
        end
    end 
     
     
    function [f, g] = costFunfmincon(v)
        Y = reshape(v, [N, rankY]);
        f = -trace(Y.'*X*Y)/2;
        if nargout > 1
            g = -X.'*Y/2 - X*Y/2;
            g = g(:);
        end
    end

    function [c, ceq, gradc, gradceq] = nonlcon(v)
        ceq = zeros(1,1);
        ceq(1,1) = v.' * v - 1;
        c = [];
        if nargout >2
            gradc = [];
            gradceq = 2 * v;
        end
    end

    function val = costfun(Y)
        val = -trace(Y.'*X*Y)/2;
    end

    function val = gradfun(Y)
        val = -X.'*Y/2 - X*Y/2;
    end


    % added the following
    function val = hessfun(Y, U)
        val = -X.'*U/2 - X*U/2;
    end

    function manvio = manifoldViolation(x)
        %Sphere Factory:
        y = x(:);
        manvio = abs(y.'*y - 1);
    end

end

