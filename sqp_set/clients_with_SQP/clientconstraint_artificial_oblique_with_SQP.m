function data = clientconstraint_artificial_oblique_with_SQP(L, rankY, methodoptions, specifier, setting)

data = NaN(3,7);
[N, ~] = size(L);
manifold = obliquefactory(rankY, N);
problem.M = manifold;
problem.cost = @(u) costFun(u);
problem.egrad = @(u) gradFun(u);
problem.ehess = @(u, d) hessFun(u, d);
x0 = problem.M.rand();

%     DEBUG only
%     checkgradient(problem);  % OK
%     checkhessian(problem);  % OK

%-------------------------Set-up Constraints-----------------------
colones = ones(N, 1);

eq_constraints_cost = cell(1,1);
eq_constraints_cost{1} = @(U) norm(U*colones,2)^2;

eq_constraints_grad = cell(1,1);
eq_constraints_grad{1} = @(U) meanzero_eq_constraint(U);

eq_constraints_hess = cell(1,1);
eq_constraints_hess{1} = @(U, D) hess_meanzero_eq_constraint(U, D);

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;
problem.eq_constraint_hess = eq_constraints_hess;

A = rand(N) * rand(N).';
A = 0.5 * (A+A.');
ineq_constraints_cost = cell(1,1);
ineq_constraints_cost{1} = @(U) trace((U) * A * (U.')) - rand(1);
ineq_constraints_grad = cell(1,1);
ineq_constraints_grad{1} = @(U) U*A + U*A.';
ineq_constraints_hess = cell(1,1);
ineq_constraints_hess{1} = @(U, D) D*A + D*A.';

problem.ineq_constraint_cost = ineq_constraints_cost;
problem.ineq_constraint_grad = ineq_constraints_grad;
problem.ineq_constraint_hess = ineq_constraints_hess;



%     Debug Only
%     checkconstraints_upto2ndorder(problem)

condet = constraintsdetail(problem);

%     ------------------------- Solving ---------------------------
    options = methodoptions;
    
    if specifier.ind(1)
        %MINI-SUM-MAX
        fprintf('Starting Mini-sum-max \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaMinimax(problem, x0, options);
        time = toc(timetic);
        
        filename = sprintf('AO_Mini-Sum-Max_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
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
        filename = sprintf('AO_ALM_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        info = rmfield(info, 'lambdas');
        info = rmfield(info, 'gammas');
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
        
        filename = sprintf('AO_LQH_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
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

        filename = sprintf('AO_LSE_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 4) = maxviolation;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
    if specifier.ind(5)
        %FMINCON interior point
        fprintf('Starting fmincon_interior_point \n');
        maxiter = methodoptions.maxOuterIter;
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon','Algorithm','interior-point', 'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
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
        
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], [], [], @nonlcon, options);
        time = toc(timetic);
        
        history.iter(1,:) =[];
        history.cost(1,:) = [];
        filename = sprintf('AO_fmincon_interior_point_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        struct2csv(info, filename);      
        
        xfinal = reshape(xfinal, [rankY, N]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 5) = output.constrviolation;
        data(2, 5) = cost;
        data(3, 5) = time;
    end
    if specifier.ind(6)
        %FMINCON sequential quadratic programming
        fprintf('Starting fmincon_SQP \n');
        maxiter = methodoptions.maxOuterIter;
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
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
        
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], [], [], @nonlcon, options);
        time = toc(timetic);
        
        history.iter(1,:) =[];
        history.cost(1,:) = [];
        filename = sprintf('AO_fmincon_SQP_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        struct2csv(info, filename);      
        
        xfinal = reshape(xfinal, [rankY, N]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 6) = output.constrviolation;
        data(2, 6) = cost;
        data(3, 6) = time;
    end     
    
    if specifier.ind(7)
        % Riemannian SQP
        fprintf('Starting Riemannian SQP \n');
        sqpoptions.maxtime = methodoptions.maxtime;
        sqpoptions.maxiter = methodoptions.maxOuterIter;
        sqpoptions.tolqpnorm = methodoptions.minstepsize;
        sqpoptions.tolgradnorm = methodoptions.minstepsize;
        sqpoptions.toliterdist = methodoptions.minstepsize;
        sqpoptions.verbosity = methodoptions.verbosity;
        timetic = tic();
        [xfinal, costfinal, info, options] = SQP(problem, x0, sqpoptions);
        time = toc(timetic);
        filename = sprintf('AO_Riemannian_SQP_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        struct2csv(info, filename);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 7) = maxviolation;
        data(2, 7) = cost;
        data(3, 7) = time;
    end
    
        filename = sprintf('AO_Info_nrep%dDim%dDen%.3f.csv',setting.repeat,setting.dim, setting.density);
        struct2csv(setting, filename);
    
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
        Y = reshape(v, [rankY, N]);
        f = trace((Y) * L * (Y.') );
        if nargout > 1
            g = Y*L + Y*L.';
            g = g(:);
        end
    end 

    function [c, ceq, gradc, gradceq] = nonlcon(v)
        Y = reshape(v, [rankY, N]);
        ceq = zeros(N+1,1);
        for rowCeq = 1: N
            ceq(rowCeq,1) = Y(:,rowCeq).'*Y(:,rowCeq) - 1;
        end
        ceq(N+1,1) = colones.' *(Y.')* Y * colones;
        c = [];
        if nargout > 2
            gradc = [];
            gradceq = zeros(rankY*N, N+1);
            grad = 2* repmat(Y*colones, 1, N);
            gradceq(:, N+1) = grad(:);
            for rowCeq = 1: N
                grad = zeros(rankY, N);
                grad(:, rowCeq) = 2*Y(:, rowCeq);
                gradceq(:, rowCeq) = grad(:);
            end
        end
    end

    function val = meanzero_eq_constraint(U)
        val = 2* repmat(U*colones, 1, N);
    end

    function val = hess_meanzero_eq_constraint(U, D)
        val = 2* repmat(D*colones, 1, N);
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
        %Oblique Factory:
        manvio = max(abs(diag(x.'*x)-colones));
    end
end

