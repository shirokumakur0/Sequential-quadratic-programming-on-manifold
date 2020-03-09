function data = clientconstraint_sphere_nonnegativePCA_with_SQP(X, rankY, methodoptions, specifier, setting)
%Y is a square matrix
%rank is smaller than columns of Y
[N, ~] = size(X);

data = NaN(3, 6);

M = spherefactory(N, rankY);
problem.M = M;
problem.cost = @(u) costfun(u);
problem.egrad = @(u) gradfun(u);
problem.ehess = @(u, d) hessfun(u, d);
x0 = M.rand();

%     DEBUG only
%     checkgradient(problem);
%     checkhessian(problem); % almost OK

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
        constrainthess = zeros(N, rankY);
        constrainthess(row, col) = 0;
        constraints_hess{(col-1)*N + row} = @(X, U) constrainthess;
    end
end

problem.ineq_constraint_cost = constraints_cost;
problem.ineq_constraint_grad = constraints_grad;
problem.ineq_constraint_hess = constraints_hess;

%     Debug Only
%     checkconstraints_upto2ndorder(problem) % Okay!

condet = constraintsdetail(problem);

%     ------------------------- Solving ---------------------------
    options = methodoptions;
        
    if specifier.ind(1)
        %ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('NNPCA_ALM_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 1) = maxviolation;
        data(2, 1) = cost;
        data(3, 1) = time;
    end

    if specifier.ind(2)
        %LQH
        fprintf('Starting LQH \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglqh(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('NNPCA_LQH_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(info, filename);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 2) = maxviolation;
        data(2, 2) = cost;
        data(3, 2) = time;
    end
    
    
    if specifier.ind(3)
        %LSE
        fprintf('Starting LSE \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglse(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('NNPCA_LSE_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 3) = maxviolation;
        data(2, 3) = cost;
        data(3, 3) = time;
    end
    
    if specifier.ind(4)
        %FMINCON Interior point method
        fprintf('Starting fmincon_interior_point \n');
        maxiter = methodoptions.maxOuterIter;
        maxFuniter = 1e+50;
        fmincontolerance = 0;  % disabling that stopping condition.
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            fminconoptions = optimoptions('fmincon', 'Algorithm','interior-point','MaxIter', maxiter, 'MaxFunEvals', maxFuniter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', fmincontolerance, 'TolCon', fmincontolerance, 'TolFun', fmincontolerance);
        else
            % Use this otherwise
            fminconoptions = optimoptions('fmincon','Algorithm','interior-point', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxFuniter,...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
                'StepTolerance', fmincontolerance, 'ConstraintTolerance', fmincontolerance, 'OptimalityTolerance', fmincontolerance);
        end
        timetic = tic();
        history = struct();
        history.iter = [];
        history.cost = [];
        history.LagGradNorm = [];
        history.maxviolation = [];
        history.KKT_residual = [];
        
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, fminconoptions);
        time = toc(timetic);

        history.iter(1,:) =[];
        history.cost(1,:) = [];
        history.LagGradNorm(1,:) = [];
        history.maxviolation(1,:) = [];
        history.KKT_residual(1,:) = [];

        filename = sprintf('NNPCA_fmincon_interior_point_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(history, filename);              

        xfinal = reshape(xfinal, [N, rankY]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 4) = output.constrviolation;
        data(2, 4) = cost;
        data(3, 4) = time;
    end

    if specifier.ind(5)
        %FMINCON sequential qudratic programming
        fprintf('Starting fmincon_SQP \n');
        maxiter = methodoptions.maxOuterIter;
        maxFuniter = 1e+50;        
        fmincontolerance = 0;  % disabling the other stopping conditions.
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            fminconoptions = optimoptions('fmincon', 'Algorithm', 'sqp' ,'MaxIter', maxiter, 'MaxFunEvals', maxFuniter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', fmincontolerance, 'TolCon', fmincontolerance, 'TolFun', fmincontolerance);
        else
            % Use this otherwise
            fminconoptions = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxFuniter,...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
                'StepTolerance', fmincontolerance, 'ConstraintTolerance', fmincontolerance, 'OptimalityTolerance', fmincontolerance);
        end
        timetic = tic();
        history = struct();
        history.iter = [];
        history.cost = [];
        history.LagGradNorm = [];
        history.maxviolation = [];
        history.KKT_residual = [];
        
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, fminconoptions);
        time = toc(timetic);
        history.iter(1,:) =[];
        history.cost(1,:) = [];
        history.LagGradNorm(1,:) = [];
        history.maxviolation(1,:) = [];
        history.KKT_residual(1,:) = [];
        
        filename = sprintf('NNPCA_fmincon_SQP_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(history, filename);              
        xfinal = reshape(xfinal, [N, rankY]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 5) = output.constrviolation;
        data(2, 5) = cost;
        data(3, 5) = time;
    end
    
    if specifier.ind(6)
        % Riemannian SQP
        fprintf('Starting Riemannian SQP \n');
        timetic = tic();
        [xfinal, costfinal, info,~] = SQP(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('NNPCA_Riemannian_SQP_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 6) = maxviolation;
        data(2, 6) = cost;
        data(3, 6) = time;
    end
    
        filename = sprintf('NNPCA_Info_nrep%dDim%dSNR%.2fDel%.2fTol%d.csv',setting.repeat,setting.dim, setting.snr,setting.delta, setting.tolKKTres);
    struct2csv(setting, filename);
    
%------------------------sub functions-----------     
     
     function stop = outfun(x, optimValues, state)
        KKT_residual = optimValues.firstorderopt^2;

        [xCurc, xCurceq, ~, ~] = nonlcon(x);
        
        [cnum, ~] = size(xCurc);
        
        for indexc = 1:cnum
            cost_at_x = max(0, xCurc(indexc));
            KKT_residual = KKT_residual + cost_at_x^2;
        end
        
        [ceqnum, ~] = size(xCurceq);
        
        for indexceq = 1:ceqnum
            KKT_residual = KKT_residual + (xCurceq(indexceq))^2;
        end
        
        KKT_residual = sqrt(KKT_residual);
        
        stop = false;
        if toc(timetic) > methodoptions.maxtime
            fprintf("Time limit exceeded\n")
            stop = true;
        elseif KKT_residual < methodoptions.tolKKTres
            fprintf("KKT residual tolerance reached\n")
            stop = true;
        end
        
        history.iter = [history.iter; optimValues.iteration];
        history.cost = [history.cost;optimValues.fval];
        history.LagGradNorm = [history.LagGradNorm; optimValues.firstorderopt];
        history.maxviolation = [history.maxviolation; optimValues.constrviolation];
        history.KKT_residual = [history.KKT_residual; KKT_residual];
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

    function val = hessfun(Y, U)
        val = -X.'*U/2 - X*U/2;
    end

    function manvio = manifoldViolation(x)
        %Sphere Factory:
        y = x(:);
        manvio = abs(y.'*y - 1);
    end
end

