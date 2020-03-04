function data = clientconstraint_rank_constraints_lc_with_SQP(m, n, k, lc_dim, P, PA, methodoptions, specifier, setting)
%% Manifold factory

data = NaN(3, 7);

M = fixedrankembeddedfactory(m, n, k); % m: # of rows, n: # of cols, k: rank
problem.M = M;
x0 = M.rand();

% Define the problem cost function. The input X is a structure with
% fields U, S, V representing a rank k matrix as U*S*V'.
% f(X) = 1/2 * || P.*(X-A) ||^2
problem.cost = @objcost;
function f = objcost(X)
    % Note that it is very much inefficient to explicitly construct the
    % matrix X in this way. Seen as we only need to know the entries
    % of Xmat corresponding to the mask P, it would be far more
    % efficient to compute those only.
    Xmat = X.U*X.S*X.V';
    f = .5*norm( P.*Xmat - PA , 'fro')^2;
end

% Define the Euclidean gradient of the cost function, that is, the
% gradient of f(X) seen as a standard function of X.
% nabla f(X) = P.*(X-A)
problem.egrad = @eobjgrad;
function G = eobjgrad(X)
    % Same comment here about Xmat.
    Xmat = X.U*X.S*X.V';
    G = P.*Xmat - PA;
end

% This is optional, but it's nice if you have it.
% Define the Euclidean Hessian of the cost at X, along H, where H is
% represented as a tangent vector: a structure with fields Up, Vp, M.
% This is the directional derivative of nabla f(X) at X along Xdot:
% nabla^2 f(X)[Xdot] = P.*Xdot
problem.ehess = @euclidean_objhessian;
function ehess = euclidean_objhessian(X, H)
    % The function tangent2ambient transforms H (a tangent vector) into
    % its equivalent ambient vector representation. The output is a
    % structure with fields U, S, V such that U*S*V' is an mxn matrix
    % corresponding to the tangent vector H. Note that there are no
    % additional guarantees about U, S and V. In particular, U and V
    % are not orthonormal.
    ambient_H = problem.M.tangent2ambient(X, H);
    Xdot = ambient_H.U*ambient_H.S*ambient_H.V';
    % Same comment here about explicitly constructing the ambient
    % vector as an mxn matrix Xdot: we only need its entries
    % corresponding to the mask P, and this could be computed
    % efficiently.
    ehess = P.*Xdot;
end

%     DEBUG only
%     figure;
%     checkgradient(problem);
%     figure;
%     checkhessian(problem); % OK?(output is 4 when it should be 3)

%% Set-up Constraints

%%% Inequality constraints
%%%% Non-negativity (ineq)
%nn_constraints_cost = cell(m * n, 1);
%for row = 1: m
%    for col = 1: n
%        nn_constraints_cost{(col-1)*m + row} = @(Y) nncostfun(Y, row, col);
%    end
%end

%nn_constraints_grad = cell(m * n, 1);
%for row = 1:m
%    for col = 1:n
%        constraintgrad = zeros(m, n);
%        constraintgrad(row, col) = -1;
%        nn_constraints_grad{(col-1)*m + row} = @(U) constraintgrad;
%    end
%end

%nn_constraints_hess = cell(m * n, 1);
%for row = 1: m
%    for col = 1: n
%        constrainthess = zeros(m, n);
%        nn_constraints_hess{(col-1)*m + row} = @(X, U) constrainthess;
%    end
%end

%function val = nncostfun(Y, row, col)
%    F = Y.U * Y.S * (Y.V).';
%    val = -F(row, col);
%end

%%%% Quadratic form (ineq)
%Amat = randi([1, 5],n);
%Amat = 0.5 * (Amat + Amat');
%Amat = 0.5 * (Amat + Amat');
%setting.A = Amat;
%r = rand(1) * 10;
%setting.r = r;
%qf_constraints_cost = cell(1,1);
%qf_constraints_cost{1} = @(Y) qfcostfun(Y);
%qf_constraints_grad = cell(1,1);
%qf_constraints_grad{1} = @(Y) qfgradfun(Y);
%qf_constraints_hess = cell(1,1);
%qf_constraints_hess{1} = @(U, D) qfhessfun(U, D);

%function val = qfcostfun(Y)
%   F = Y.U * Y.S * (Y.V).';
%   val = trace((F) * Amat * (F.')) - r;
%end

%function egrad = qfgradfun(Y)
%    F = Y.U * Y.S * (Y.V).';
%    egrad = F*Amat + F*Amat.';
%end

%function ehess = qfhessfun(Y,D)
%    ambD = problem.M.tangent2ambient(Y, D);
%    F = ambD.U * ambD.S * (ambD.V).';
%    ehess = F*Amat + F*Amat.';
%end

%%% Equality constraints
%%%% Linear constraints (eq)
while true
    C = randi([1, 5], m*n, lc_dim);
    C = orth(C);
    [~,cC] = size(C);
    if cC == lc_dim
        break;
    end
end
C = C';
d = rand(lc_dim , 1) * 2;

lc_constraints_cost = cell(lc_dim, 1);
lc_constraints_grad = cell(lc_dim, 1);
lc_constraints_hess = cell(lc_dim, 1);

for row = 1 : lc_dim
    rowC = C(row,:);
    lc_constraints_cost{row} = @(Y) lccostfun(Y, rowC, row);
    reshapedrowC = reshape(rowC', [m,n]);
    lc_constraints_grad{row} = @(U) reshapedrowC;
    lc_constraints_hess{row} = @(U, D) zeros(m, n);
end

function val = lccostfun(Y, rowC, row)
    F = Y.U * Y.S * (Y.V).';
    val = rowC* F(:) -d(row);
end

% condet constraints
%ineq_constraints_cost = [nn_constraints_cost; qf_constraints_cost];
%ineq_constraints_grad = [nn_constraints_grad; qf_constraints_grad];
%ineq_constraints_hess = [nn_constraints_hess; qf_constraints_hess];

eq_constraints_cost = lc_constraints_cost;
eq_constraints_grad = lc_constraints_grad;
eq_constraints_hess = lc_constraints_hess;

% constraints setting
%problem.ineq_constraint_cost = ineq_constraints_cost;
%problem.ineq_constraint_grad = ineq_constraints_grad;
%problem.ineq_constraint_hess = ineq_constraints_hess;

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;
problem.eq_constraint_hess = eq_constraints_hess;


%     Debug Only
%     checkconstraints_upto2ndorder(problem) % OK...? (qfocnstraints:  output is 4 when it should be 3)


condet = constraintsdetail(problem);

%% Calculating by solvers

    options = methodoptions;
        
    % The solver is not available.
    %if specifier.ind(1)
    %    %MINI-SUM-MAX
    %    fprintf('Starting Mini-sum-max \n');
    %    timetic = tic();
    %    [xfinal, info] = fixedrank_exactpenaltyViaMinimax(problem, x0, options);
    %    time = toc(timetic);
    %    filename = sprintf('RC_Mini-Sum-Max_nrep%dRowdim%dColdim%dRank%dDegnum%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.deg_num);
    %    struct2csv(info, filename);

    %    [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
    %    % maxviolation = max(maxviolation, manifoldViolation(xfinal));
    %    data(1, 1) = maxviolation;
    %    data(2, 1) = cost;
    %    data(3, 1) = time;
    %end
    
    if specifier.ind(2)
        %ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_lc_ALM_nrep%dRowdim%dColdim%dRank%dLcdim%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.lc_dim);
        % info = rmfield(info,'lambdas');
        % info = rmfield(info,'gammas');
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        % maxviolation = max(maxviolation, manifoldViolation(xfinal));
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
        filename = sprintf('RC_lc_LQH_nrep%dRowdim%dColdim%dRank%dLcdim%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.lc_dim);
        struct2csv(info, filename);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        % maxviolation = max(maxviolation, manifoldViolation(xfinal));
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
        filename = sprintf('RC_lc_LSE_nrep%dRowdim%dColdim%dRank%dLcdim%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.lc_dim);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        % maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 4) = maxviolation;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
 %   if specifier.ind(5)
 %       %FMINCON Interior point method
 %       fprintf('Starting fmincon_interior_point \n');
 %       maxiter = methodoptions.maxOuterIter;    
 %       if specifier.matlabversion == 0
 %           % Use this if you are at 2015a or older.
 %           options = optimoptions('fmincon', 'Algorithm','interior-point','MaxIter', maxiter, 'MaxFunEvals', maxiter,...
 %               'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
 %               'TolX', methodoptions.minstepsize);
 %       else
 %           % Use this otherwise
 %           options = optimoptions('fmincon','Algorithm','interior-point', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
 %               'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
 %               'StepTolerance', methodoptions.minstepsize);
 %       end
 %       timetic = tic();
 %       history = struct();
 %       history.iter = [];
 %       history.cost = [];
 %       %history.maxviolation = [];
 %       %history.meanviolation = [];
 %       
 %       [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, options);
 %       time = toc(timetic);
 %
 %       history.iter(1,:) =[];
 %       history.cost(1,:) = [];
 %       filename = sprintf('RC_fmincon_interior_point_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
 %       struct2csv(history, filename);
 %
 %       xfinal = reshape(xfinal, [N, rankY]);
 %       [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
 %       data(1, 5) = output.constrviolation;
 %       data(2, 5) = cost;
 %       data(3, 5) = time;
 %   end

 %   if specifier.ind(6)
 %       %FMINCON sequential qudratic programming
 %       fprintf('Starting fmincon_interior_point \n');
 %       maxiter = methodoptions.maxOuterIter;    
 %       if specifier.matlabversion == 0
 %           % Use this if you are at 2015a or older.
 %           options = optimoptions('fmincon', 'Algorithm', 'sqp' ,'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
 %               'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
 %               'TolX', methodoptions.minstepsize);
 %       else
 %           % Use this otherwise
 %           options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
 %               'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
 %               'StepTolerance', methodoptions.minstepsize);
 %       end
 %       timetic = tic();
 %       history = struct();
 %       history.iter = [];
 %       history.cost = [];
 %       % history.maxviolation = [];
 %       % history.meanviolation = [];
 %       [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], zeros(N*rankY, 1), [], @nonlcon, options);
 %       time = toc(timetic);
 %       history.iter(1,:) =[];
 %       history.cost(1,:) = [];
 %       filename = sprintf('RC_fmincon_SQP_nrep%dDim%dSNR%.2fDel%.2f.csv',setting.repeat,setting.dim, setting.snr,setting.delta);
 %       struct2csv(history, filename);              
 %       xfinal = reshape(xfinal, [N, rankY]);
 %       [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
 %       data(1, 6) = output.constrviolation;
 %       data(2, 6) = cost;
 %       data(3, 6) = time;
 %   end
    
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
        filename = sprintf('RC_lc_Riemannian_SQP_nrep%dRowdim%dColdim%dRank%dLcdim%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.lc_dim);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        % maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1, 7) = maxviolation;
        data(2, 7) = cost;
        data(3, 7) = time;
    end
    
    filename = sprintf('RC_lc_Info_nrep%dRowdim%dColdim%dRank%dLcdim%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.lc_dim);
    struct2csv(setting, filename);
    
%% Sub functions
%%% for fmincon (not maintained yet)
     function stop = outfun(x, optimValues, state)
        stop = false;
        if toc(timetic) > methodoptions.maxtime
            stop = true;
        end
        history.iter = [history.iter; optimValues.iteration];
        history.cost = [history.cost;optimValues.fval];
        % [history.maxviolation, history.meanviolation, cost] = evaluation(problem, x, condet);
    end 
     
     
    function [f, g] = costFunfmincon(v)
        Y = reshape(v, [N, rankY]);
        f = -trace(Y.' *X * Y)/2;
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

%    function val = costfun(Y)
%        F = Y.U * Y.S * (Y.V).';
%        val = -0.5 * trace(F.' * X * F);
%    end
%
%    function val = gradfun(Y)
%        F = Y.U * Y.S * (Y.V).';
%        val = -0.5 * X.'*F - 0.5 *X*F;
%    end
%
%    function val = hessfun(Y, D)
%        % The function tangent2ambient transforms H (a tangent vector) into
%        % its equivalent ambient vector representation. The output is a
%        % structure with fields U, S, V such that U*S*V' is an mxn matrix
%        % corresponding to the tangent vector H. Note that there are no
%        % additional guarantees about U, S and V. In particular, U and V
%        % are not orthonormal.
%        ambD = problem.M.tangent2ambient(Y, D);
%        F = ambD.U * ambD.S * (ambD.V).';
%        val = - 0.5 * X.'*F - 0.5 * X*F;
%    end

    % function manvio = manifoldViolation(x)
    %     y = x(:);
    %     manvio = abs(y.'*y - 1);
    % end
end