function data = clientconstraint_fixedrank_equality_nonneg_lowrankcompletion(m, n, k, P, A, eqindices, methodoptions, specifier, setting)
%% Manifold factory

data = NaN(3, 4);

M = fixedrankembeddedfactory(m, n, k); % m: # of rows, n: # of cols, k: rank
problem.M = M;

% initial point will be created later

%% Setting the file path
filepath = setting.filepath;
%% Set-up Constraints

%%% nonnegativity (ineq) and equality (eq)
[eqnum,~] = size(eqindices);
ineqnum = m*n - eqnum;

eq_constraints_cost = cell(eqnum,1);
eq_constraints_grad = cell(eqnum,1);
eq_constraints_hess = cell(eqnum,1);
nn_constraints_cost = cell(ineqnum, 1);
nn_constraints_grad = cell(ineqnum, 1);
nn_constraints_hess = cell(ineqnum, 1);
eqconst_idx = 1;
nnconst_idx = 1;

for row = 1:m
    for col = 1:n
        num = (col - 1) * m + row;
        flag = find(eqindices == num,1);
        if ~isempty(flag)
            eq_constraints_cost{eqconst_idx} = @(Y) eqcostfun(Y, row, col);
            
            constraintgrad = zeros(m, n);
            constraintgrad(row, col) = 1;
            eq_constraints_grad{eqconst_idx} = @(U) constraintgrad;
            
            constrainthess = zeros(m, n);
            eq_constraints_hess{eqconst_idx} = @(X, U) constrainthess;
            eqconst_idx = eqconst_idx + 1;
        else
            nn_constraints_cost{nnconst_idx} = @(Y) nncostfun(Y, row, col);
            
            constraintgrad = zeros(m, n);
            constraintgrad(row, col) = -1;
            nn_constraints_grad{nnconst_idx} = @(U) constraintgrad;
            
            constrainthess = zeros(m, n);
            nn_constraints_hess{nnconst_idx} = @(X, U) constrainthess;
            nnconst_idx = nnconst_idx + 1;
        end
    end
end

function val = nncostfun(Y, row, col)
    Vt = Y.V.';
    val = - Y.U(row,:) * Y.S * Vt(:,col);
end

function val = eqcostfun(Y, row, col)
    Vt = Y.V.';
    val = Y.U(row,:) * Y.S * Vt(:,col) - A(row,col);
end

ineq_constraints_cost = nn_constraints_cost;
ineq_constraints_grad = nn_constraints_grad;
ineq_constraints_hess = nn_constraints_hess;

% constraints setting
problem.ineq_constraint_cost = ineq_constraints_cost;
problem.ineq_constraint_grad = ineq_constraints_grad;
problem.ineq_constraint_hess = ineq_constraints_hess;

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;
problem.eq_constraint_hess = eq_constraints_hess;

%     Debug Only
%     checkconstraints_upto2ndorder(problem) 
     
condet = constraintsdetail(problem);

%% Setting a objective function
% Note that the observed elements (used as equality constraints) are
% eliminated when constructing the objective function.

% Define the problem cost function. The input X is a structure with
% fields U, S, V representing a rank k matrix as U*S*V'.
% f(X) = 1/2 * || P.*(X-A) ||^2
PA = P .* A;

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
%      figure;
%      checkgradient(problem);
%      figure;
%      checkhessian(problem); 


%% Generating x0
if strcmp(setting.initialpoint, "eye")
    x0 = struct();
    x0.U = [eye(k);zeros(m-k,k)];
    x0.S = eye(k);
    x0.V = [eye(k);zeros(n-k,k)];
    
elseif strcmp(setting.initialpoint, "feasible_region")
    feasible_problem.M = M;
    feasible_problem.cost = @zerofun;
    feasible_problem.egrad = @egradzerofun;
    feasible_problem.ehess = @ehesszerofun;
    feasible_problem.ineq_constraint_cost = ineq_constraints_cost;
    feasible_problem.ineq_constraint_grad = ineq_constraints_grad;
    feasible_problem.ineq_constraint_hess = ineq_constraints_hess;
    feasible_problem.eq_constraint_cost = eq_constraints_cost;
    feasible_problem.eq_constraint_grad = eq_constraints_grad;
    feasible_problem.eq_constraint_hess = eq_constraints_hess; 
    feasible_x0 = M.rand();
    feasible_options.maxOuterIter = 200;
    feasible_options.maxtime = 40;  % 40
    feasible_options.outerverbosity = 1;  % 1
    feasible_options.tolKKTres = 10^(-2);  
    feasible_options.startingtolgradnorm = 1;
    feasible_options.endingtolgradnorm = feasible_options.tolKKTres;
    
    %     Debug Only
    %     checkconstraints_upto2ndorder(problem) 
    
    fprintf('Starting LQH to calculate a feasible point\n');
    [x0, info, ~] = exactpenaltyViaSmoothinglqh(feasible_problem, feasible_x0, feasible_options);
    filename = sprintf('RC_nnlc_LQH_feasible_initial_point_%s.csv',filepath);
    struct2csv(info, filename);      
else
    x0 = M.rand();
end
setting.x0 = x0.U * x0.S * x0.V';

%% Calculating by solvers

    options = methodoptions;
    

    if specifier.ind(1)
        % ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info, residual] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_nnlc_ALM_%s.csv',filepath);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        % DEBUG only
        % checkSVD(xfinal);
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
        filename = sprintf('RC_nnlc_LQH_%s.csv',filepath);
        struct2csv(info, filename);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        % DEBUG only
        % checkSVD(xfinal);
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
        filename = sprintf('RC_nnlc_LSE_%s.csv',filepath);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        % DEBUG only
        % checkSVD(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        data(1, 3) = residual;
        data(2, 3) = cost;
        data(3, 3) = time;
    end
    
    if specifier.ind(4)
        % Riemannian SQP
        fprintf('Starting Riemannian SQP \n');
        timetic = tic();
        [xfinal, costfinal, residual, info,~] = SQP(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_nnlc_Riemannian_SQP_%s.csv',filepath);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        % DEBUG only
        % checkSVD(xfinal);
        data(1, 4) = residual;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
    filename = sprintf('RC_nnlc_Info_%s.csv',filepath);
    struct2csv(setting, filename);
    
%% Sub functions
    function rankflag = check_rank(xfinal)
        xfinalmat = xfinal.U * xfinal.S * xfinal.V';
        r = rank(xfinalmat);
        rankflag = 0;
        if r ~= k
            fprintf("The rank of the output is %d whereas required one is %d\n", r, k) 
        else
            fprintf("Rank constraint satisfied\n")
            rankflag = 1;
        end
    end

    function checkSVD(xfinal)
       xfinalmat = xfinal.U * xfinal.S * xfinal.V';
       e = svd(xfinalmat)
    end

    function ans = zerofun(X)
        ans = 0;
    end

    function ans = egradzerofun(X)
        ans = zeros(m,n);
    end

    function ans = ehesszerofun(X,H)
        ans = zeros(m,n);
    end
end