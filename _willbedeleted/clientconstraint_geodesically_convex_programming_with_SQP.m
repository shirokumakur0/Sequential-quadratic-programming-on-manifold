function data = clientconstraint_geodesically_convex_programming_with_SQP(n, m, methodoptions, specifier, setting)
%n: dimension
%m: size (the number of matrices)
%K: max iterations
[N, ~] = size(X);

data = NaN(3, 7);

%% Generate collection of PSD
A=genPosdef(n,m);  % A: cells data, whose number of m, each of which has a n*n sym. mat. 
cond_no = cond_numbers(A);

% generate weights
w=rand(1,m);
s=sum(w);
w=w./s;

% means
am=arithmeticMean(A);
[hm Ai] =harmonicMean(A);

M = sympositivedefinitefactory_mod(n);
problem.M = M;

% initialization
means{1}=am;
means{2}=hm;
x0 = hm;

%% already imported the cadidate functions of cost grad and hess?
 % firstly, check the validity of them!
 
problem.cost = @(u) costfun(u);
problem.egrad = @(u) gradfun(u);
problem.ehess = @(u, d) hessfun(u, d);
x0 = M.rand();

%     DEBUG only
     checkgradient(problem);
     checkhessian(problem); 

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
        %MINI-SUM-MAX
        fprintf('Starting Mini-sum-max \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaMinimax(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('NNPCA_Mini-Sum-Max_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
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
        filename = sprintf('NNPCA_ALM_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
        % info = rmfield(info,'lambdas');
        % info = rmfield(info,'gammas');
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
        filename = sprintf('NNPCA_LQH_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
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
        filename = sprintf('NNPCA_LSE_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 4) = maxviolation;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
    if specifier.ind(5)
        %FMINCON Interior point method
        fprintf('Starting fmincon_interior_point \n');
        maxiter = methodoptions.maxOuterIter;    
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon', 'Algorithm','interior-point','MaxIter', maxiter, 'MaxFunEvals', maxiter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', methodoptions.minstepsize);
        else
            % Use this otherwise
            options = optimoptions('fmincon','Algorithm','interior-point', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
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

        history.iter(1,:) =[];
        history.cost(1,:) = [];
        filename = sprintf('NNPCA_fmincon_interior_point_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
        struct2csv(history, filename);              

        xfinal = reshape(xfinal, [N, rankY]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1, 5) = output.constrviolation;
        data(2, 5) = cost;
        data(3, 5) = time;
    end

    if specifier.ind(6)
        %FMINCON sequential qudratic programming
        fprintf('Starting fmincon_interior_point \n');
        maxiter = methodoptions.maxOuterIter;    
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon', 'Algorithm', 'sqp' ,'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
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
        history.iter(1,:) =[];
        history.cost(1,:) = [];
        filename = sprintf('NNPCA_fmincon_SQP_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
        struct2csv(history, filename);              
        xfinal = reshape(xfinal, [N, rankY]);
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
        sqpoptions.mineigval_correction = methodoptions.mineigval_correction;

        timetic = tic();
        [xfinal, costfinal, info,~] = SQP(problem, x0, sqpoptions);
        time = toc(timetic);
        filename = sprintf('NNPCA_Riemannian_SQP_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 7) = maxviolation;
        data(2, 7) = cost;
        data(3, 7) = time;
    end
    
    filename = sprintf('NNPCA_Info_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
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

    function f = costfun(x, A)  % cost function
      if ~isposdef(A)
        f = Inf;
      else
        f = gmobj(x,A);
      end
    end

    function g = egradfun(x, A) % gradient of f(x)
      if ~isposdef(A)
          g  = Inf(size(x));
          %gf = Inf(size(x));
      else
          g = zeros(size(x));
          for i=1:numel(A)
            g = g + logm(A{i}\x);
          end
          g=2*g/x;
          g=(g+g')/2;
          %gf = problem.M.egrad2rgrad(x, g);
      end
    end

    function output = hessfun(x, v, A) % TODO CONSIDER THE VALIDITY
    k = length(A.inv);
    output.TV = zeros(size(x.U));
    fv = full_pd(x, v);
    for i = 1:k
        temp = fv.TV * A.invh{i} * x.log{i} * A.h{i};
        output.TV = output.TV + temp - temp' + ...
            x.U * Dlog(A.inv{i} * x.U, A.inv{i} * fv.TV);
    end
    output.TV = 0.5 * (output.TV + output.TV');
    output.TV = output.TV/k;
    output = vech_pd(x, output);
    end

    %function val = costfun(Y)
    %    val = -trace(Y.'*X*Y)/2;
    %end

    %function val = gradfun(Y)
    %    val = -X.'*Y/2 - X*Y/2;
    %end

    %function val = hessfun(Y, U)
    %    val = -X.'*U/2 - X*U/2;
    %end


    function manvio = manifoldViolation(x)
        %Sphere Factory:
        y = x(:);
        manvio = abs(y.'*y - 1);
    end
end
