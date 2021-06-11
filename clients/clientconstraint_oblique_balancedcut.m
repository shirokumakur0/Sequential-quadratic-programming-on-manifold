function data = clientconstraint_oblique_balancedcut(L, rankY, methodoptions, specifier, setting)

data = NaN(3,6);
[N, ~] = size(L);
manifold = obliquefactory(rankY, N);
problem.M = manifold;
problem.cost = @(u) costFun(u);
problem.egrad = @(u) gradFun(u);
problem.ehess = @(u, d) hessFun(u, d);
x0 = problem.M.rand();
setting.x0 = x0;

%     DEBUG only
%     checkgradient(problem);
%     checkhessian(problem);

%-------------------------Set-up Constraints-----------------------
colones = ones(N, 1);
[~,xC] = size(x0);
eq_gradmat1 = [ones(1,xC); zeros(1,xC)];
eq_gradmat2 = [zeros(1,xC); ones(1,xC)];

eq_constraints_cost = cell(2,1);
eq_constraints_cost{1} = @(U) U(1,:) * colones;
eq_constraints_cost{2} = @(U) U(2,:) * colones;

eq_constraints_grad = cell(2,1);
eq_constraints_grad{1} = @(U) eq_gradmat1;
eq_constraints_grad{2} = @(U) eq_gradmat2;

eq_constraints_hess = cell(2,1);
eq_constraints_hess{1} = @(U, D) 0;
eq_constraints_hess{2} = @(U, D) 0;

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;
problem.eq_constraint_hess = eq_constraints_hess;

%     Debug Only
%     checkconstraints_upto2ndorder(problem)

condet = constraintsdetail(problem);

%     ------------------------- Solving ---------------------------
    options = methodoptions;
    
    if specifier.ind(1)
        % ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info, residual] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('BC_ALM_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
        struct2csv(info, filename);

        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        
        data(1, 1) = residual;
        data(2, 1) = cost;
        data(3, 1) = time;
    end

    if specifier.ind(2)
        % LQH
        fprintf('Starting LQH \n');
        timetic = tic();
        [xfinal, info, residual] = exactpenaltyViaSmoothinglqh(problem, x0, options);
        time = toc(timetic);
        
        filename = sprintf('BC_LQH_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 2) = residual;
        data(2, 2) = cost;
        data(3, 2) = time;
    end
    
    
    if specifier.ind(3)
        % LSE
        fprintf('Starting LSE \n');
        timetic = tic();
        [xfinal, info, residual] = exactpenaltyViaSmoothinglse(problem, x0, options);
        time = toc(timetic);

        filename = sprintf('BC_LSE_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 3) = residual;
        data(2, 3) = cost;
        data(3, 3) = time;
    end
    
    if specifier.ind(4)
        % FMINCON sequential quadratic programming
        fprintf('Starting fmincon_SQP \n');
        maxiter = methodoptions.maxOuterIter;
        maxFuniter = 1e+10;        
        fmincontolerance = 0;  % disabling the other stopping conditions.        
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            fminconoptions = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', maxiter, 'MaxFunEvals', maxFuniter,...
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
        history.time = [];
        history.LagGradNorm = [];
        history.maxviolation = [];
        history.KKT_residual = [];
        
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], [], [], @nonlcon, fminconoptions);
        time = toc(timetic);
        
        history.iter(1,:) =[];
        history.cost(1,:) = [];
        history.time(1,:) = [];
        history.LagGradNorm(1,:) = [];
        history.maxviolation(1,:) = [];
        history.KKT_residual(1,:) = [];
        
        filename = sprintf('BC_fmincon_SQP_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
        struct2csv(history, filename);      
        
        xfinal = reshape(xfinal, [rankY, N]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        % data(1, 4) is blank
        data(2, 4) = cost;
        data(3, 4) = time;
    end     
    
    if specifier.ind(5)
        % Riemannian SQP
        fprintf('Starting Riemannian SQP \n');
        timetic = tic();
        [xfinal, costfinal, residual, info, ~] = SQP(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('BC_Riemannian_SQP_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 5) = residual;
        data(2, 5) = cost;
        data(3, 5) = time;
    end
    
        filename = sprintf('BC_Info_nrep%dDim%dDen%.3fTol%d.csv',setting.repeat,setting.dim, setting.density, setting.tolKKTres);
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
        elseif KKT_residual <= methodoptions.tolKKTres
            fprintf("KKT residual tolerance reached\n")
            stop = true;
        end
        
        fmincontime = toc(timetic);
        
        history.iter = [history.iter; optimValues.iteration];
        history.cost = [history.cost;optimValues.fval];
        history.time = [history.time; fmincontime];
        history.LagGradNorm = [history.LagGradNorm; optimValues.firstorderopt];
        history.maxviolation = [history.maxviolation; optimValues.constrviolation];
        history.KKT_residual = [history.KKT_residual; KKT_residual];
    end 
     
    function [f, g] = costFunfmincon(v)
        Y = reshape(v, [rankY, N]);
        f = trace((Y) * L * (Y.') );
        if nargout > 1
            g = Y*L + Y*L.';
            g = g(:);
        end
    end 

    function [c, ceq, gradc, gradceq] = nonlcon(v)
        Y = reshape(v, [rankY, N]);
        ceq = zeros(N+rankY,1);
        for rowCeq = 1: N
            ceq(rowCeq,1) = Y(:,rowCeq).'*Y(:,rowCeq) - 1;
        end
        for rowCeq = 1:rankY
            ceq(N+rowCeq,1) = Y(rowCeq, :) * colones;
        end
        c = [];
        if nargout > 2
            gradc = [];
            gradceq = zeros(rankY*N, N+rankY);
            for rowCeq = 1:rankY
                grad = zeros(rankY,N);
                grad(rowCeq, :) = 1;
                gradceq(:, N+rowCeq) = grad(:);
            end
            for rowCeq = 1: N
                grad = zeros(rankY, N);
                grad(:, rowCeq) = 2*Y(:, rowCeq);
                gradceq(:, rowCeq) = grad(:);
            end
        end
    end

    function val = costFun(u)
        val = trace((u) * L * (u.') );
    end

    function val = gradFun(u)
        val = u*L + u*L.';
    end

    function val = hessFun(u, d)
        val = d*L + d*L.';
    end

    function manvio = manifoldViolation(x)
        % Oblique Factory:
        manvio = max(abs(diag(x.'*x)-colones));
    end
end

