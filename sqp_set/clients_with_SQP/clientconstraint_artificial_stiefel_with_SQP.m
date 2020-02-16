 function data = clientconstraint_artificial_stiefel_with_SQP(D, rankY, methodoptions, specifier, setting)
% D has to be symmetric. 
% rankY is positive integer
% returned table:
%                Col 1: Mini-sum-max, Col 2, ALM, Col 3, lqh, Col 4, lse, Col 5, fmincon
% Row1: Maxviolation
% Row2: Cost
% Row3: Time

    data = NaN(3, 7);

    [N, ~] = size(D);
    

    M = stiefelfactory(N, rankY);
    problem.M = M;
    problem.cost = @(Y) costFun(Y,D);
    problem.egrad = @(Y) gradFun(Y,D);
    problem.ehess = @(Y,U) hessFun(Y,U,D);
    x0 = M.rand();

%     DEBUG only
%     figure;
%     checkgradient(problem); % Okay!
%     figure;
%     checkhessian(problem); % Okay!

%-------------------------Set-up Constraints-----------------------
% Nonnegativity of all entries
    ineq_constraints_cost = cell(1, N * rankY);
    for row = 1: N
        for col = 1: rankY
            ineq_constraints_cost{(col-1)*N + row} = @(Y) -Y(row, col);
        end
    end

    ineq_constraints_grad = cell(1, N * rankY);
    for row = 1: N
        for col = 1: rankY
            constraintgrad = zeros(N, rankY);
            constraintgrad(row, col) = -1;
            ineq_constraints_grad{(col-1)*N + row} = @(Y) constraintgrad;
        end
    end
    

    ineq_constraints_hess = cell(1, N * rankY);
    for row = 1: N
        for col = 1: rankY
            constrainthess = zeros(N, rankY);
            constrainthess(row, col) = 0;
            ineq_constraints_hess{(col-1)*N + row} = @(Y, U) constrainthess;
        end
    end
    
    problem.ineq_constraint_cost = ineq_constraints_cost;
    problem.ineq_constraint_grad = ineq_constraints_grad;
    problem.ineq_constraint_hess = ineq_constraints_hess;

% 1 is the eigenvector of Y*Y.'
    colones = ones(N,1);
    eq_constraints_cost = cell(1, N);
    for row = 1: N
        eq_constraints_cost{row} = @(Y) Y(row,:)*Y'*colones -1;
    end

    eq_constraints_grad = cell(1, N);
    for row = 1: N
        eq_constraints_grad{row} = @(Y) constrainteqgradFun(Y, row);
    end
    
    eq_constraints_hess = cell(1, N);
    for row = 1: N
        eq_constraints_hess{row} = @(Y, U) constrainteqhessFun(Y,U, row);
    end
    
    
    problem.eq_constraint_cost = eq_constraints_cost;
    problem.eq_constraint_grad = eq_constraints_grad;
    problem.eq_constraint_hess = eq_constraints_hess;

    condet = constraintsdetail(problem);
    
%     Debug Only
%     checkconstraints_upto2ndorder(problem) % OKAY!
    
%     ------------------------- Solving ---------------------------
    options = methodoptions;
    
    if specifier.ind(1)
        %MINI-SUM-MAX
        fprintf('Starting Mini-sum-max \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaMinimax(problem, x0, options);
        time = toc(timetic);

        filename = sprintf('AS_Mini-Sum-Max_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 1) = maxviolation;
        data(2, 1) = cost;
        data(3, 1) = time;
    end    
    
    if specifier.ind(2)
        %ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        
        info = rmfield(info, 'lambdas');
        info = rmfield(info, 'gammas');
        
        filename = sprintf('AS_ALM_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 2) = maxviolation;
        data(2, 2) = cost;
        data(3, 2) = time;
    end
    
    if specifier.ind(3)
        %LQH
        fprintf('Starting LQH \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglqh(problem, x0, options);
        time = toc(timetic);
        
        filename = sprintf('AS_LQH_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(info, filename);

        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 3) = maxviolation;
        data(2, 3) = cost;
        data(3, 3) = time;
    end
    
    
    if specifier.ind(4)
        %LSE
        fprintf('Starting LSE \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglse(problem, x0, options);
        time = toc(timetic);

        filename = sprintf('AS_LQH_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 4) = maxviolation;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
    if specifier.ind(5)
        %FMINCON interior-point
        fprintf('Starting fmincon_interior_point \n');
        maxiter = methodoptions.maxOuterIter;
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon', 'Algorithm', 'interior-point' , 'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', methodoptions.minstepsize);
        else
            % Use this otherwise
            options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
                'StepTolerance', methodoptions.minstepsize);
        end
        timetic = tic();
        history = struct();
        history.iter = [];
        history.cost = [];
        %history.maxviolation = [];
        %history.meanviolation = [];
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, options);
        time = toc(timetic);
        
        history.iter(1,:) = [];
        history.cost(1,:)  = [];
        
        filename = sprintf('AS_fmincon_interior_point_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(history, filename);        
        
        xfinal = reshape(xfinal, [N, rankY]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 5) = output.constrviolation;
        data(2, 5) = cost;
        data(3, 5) = time;
    end

    if specifier.ind(6)
        %FMINCON sequential qudratic programming
        fprintf('Starting fmincon_SQP \n');
        maxiter = methodoptions.maxOuterIter;
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon', 'Algorithm', 'sqp' , 'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', methodoptions.minstepsize);
        else
            % Use this otherwise
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
                'StepTolerance', methodoptions.minstepsize);
        end
        timetic = tic();
        history = struct();
        history.iter = [];
        history.cost = [];
        %history.maxviolation = [];
        %history.meanviolation = [];
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, options);
        time = toc(timetic);
        
        history.iter(1,:) = [];
        history.cost(1,:)  = [];
        
        filename = sprintf('AS_fmincon_SQP_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(history, filename);        
        
        xfinal = reshape(xfinal, [N, rankY]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 6) = output.constrviolation;
        data(2, 6) = cost;
        data(3, 6) = time;
    end
    
    if specifier.ind(7)
        %Riemannian SQP
        fprintf('Starting Riemannian SQP \n');
        timetic = tic();
        [xfinal, costfinal, info, options] = SQP(problem, x0, options);
        time = toc(timetic);

        filename = sprintf('AS_Riemannian_SQP_nrep%dNum%dDim%dRank%d.csv',setting.repeat,...
            setting.data_num, setting.data_dim, setting.rank);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 7) = maxviolation;
        data(2, 7) = cost;
        data(3, 7) = time;
    end
    
     %------------------------sub functions-----------     
     
    function stop = outfun(x, optimValues, state)
        stop = false;
        if toc(timetic) > methodoptions.maxtime
            stop = true;
        end
        history.iter = [history.iter; optimValues.iteration];
        history.cost = [history.cost;optimValues.fval];
        %[history.maxviolation, history.meanviolation, cost] = evaluation(problem, x, condet);
    end 
     
    
    function [f, g] = costFunfmincon(v)
        Y = reshape(v, [N, rankY]);
        f = 0;
        for ii = 1:rankY
            f = f + Y(:,ii)'*D*Y(:,ii)/2;
        end
        if nargout > 1
            g = D * Y;
            g = g(:);
        end
    end
    
    function [c, ceq, gradc, gradceq] = nonlcon(v)
        Y = reshape(v, [N, rankY]);
        Diff1 = Y'*Y-eye(rankY);
        Diff2 = Y*Y'*colones - colones;
        ceq = [Diff1(:); Diff2];
        c = [];
        if nargout >2
            gradc = [];
            gradceq = zeros(N*rankY, rankY^2+N);
            for rowceq = 1 : rankY
                for colceq = 1 : rankY
                    grad = zeros(N, rankY);
                    grad(:, colceq) = grad(:, colceq) + Y(:, rowceq);
                    grad(:, rowceq) = grad(:, rowceq) + Y(:, colceq);
                    gradceq(:, rankY * (colceq-1) + rowceq) = grad(:);
                end
            end
            for ii = 1: N
                    val = repmat(Y(ii, :), N, 1);
                    val(ii, :) = val(ii, :) + colones'*Y;
                    gradceq(:, rankY^2 + ii) = val(:);
            end
        end
    end
    

    function f = costFun(Y,D)
        f = 0;
        for ii = 1:rankY
            f = f + Y(:,ii)'*D*Y(:,ii)/2;
        end
    end

    function val = gradFun(Y,D)
        val = D*Y;
    end

    function val = hessFun(Y,U,D)
        val = D*U;
    end

    function val = constrainteqgradFun(Y, i)
        val = repmat(Y(i,:), N, 1);
        val(i, :) = val(i,:) + colones'*Y;
    end

    function val = constrainteqhessFun(Y,U, i)
        val = repmat(U(i,:), N, 1);
        val(i, :) = val(i,:) + colones'*U;
    end
    
    function displayinfo(stats)
        figure;
        subplot(2,2,1)
        plot([stats.time], [stats.maxviolation], '.-');
        xlabel('Iter');
        ylabel('Maxviolation');
        
        subplot(2,2,2)
        plot([stats.time], [stats.meanviolation], '.-');
        xlabel('Iter');
        ylabel('Meanviolation');
        
        subplot(2,2,3)
        plot([stats.time], [stats.cost], '.-');
        xlabel('Iter');
        ylabel('cost');
    end

    function manvio = manifoldViolation(x)
        %Stefel Factory:
        manvio = max(max(abs(x.'*x-eye(rankY))));
    end
end