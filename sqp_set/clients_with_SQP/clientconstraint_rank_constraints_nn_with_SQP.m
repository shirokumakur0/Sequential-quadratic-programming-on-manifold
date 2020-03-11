function data = clientconstraint_rank_constraints_nn_with_SQP(m, n, k, P, PA, methodoptions, specifier, setting)
%% Manifold factory

data = NaN(3, 4);

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
nn_constraints_cost = cell(m * n, 1);
for row = 1: m
    for col = 1: n
        nn_constraints_cost{(col-1)*m + row} = @(Y) nncostfun(Y, row, col);
    end
end

nn_constraints_grad = cell(m * n, 1);
for row = 1:m
    for col = 1:n
        constraintgrad = zeros(m, n);
        constraintgrad(row, col) = -1;
        nn_constraints_grad{(col-1)*m + row} = @(U) constraintgrad;
    end
end

nn_constraints_hess = cell(m * n, 1);
for row = 1: m
    for col = 1: n
        constrainthess = zeros(m, n);
        nn_constraints_hess{(col-1)*m + row} = @(X, U) constrainthess;
    end
end

function val = nncostfun(Y, row, col)
    Vt = Y.V.';
    val = - Y.U(row,:) * Y.S * Vt(:,col);
    
    %F = Y.U * Y.S * (Y.V).';
    %val = -F(row, col);
end


ineq_constraints_cost = nn_constraints_cost;
ineq_constraints_grad = nn_constraints_grad;
ineq_constraints_hess = nn_constraints_hess;

% constraints setting
problem.ineq_constraint_cost = ineq_constraints_cost;
problem.ineq_constraint_grad = ineq_constraints_grad;
problem.ineq_constraint_hess = ineq_constraints_hess;

%     Debug Only
%     checkconstraints_upto2ndorder(problem) % OK...? (qfocnstraints:  output is 4 when it should be 3)


condet = constraintsdetail(problem);

%% Calculating by solvers

    options = methodoptions;
        
    if specifier.ind(1)
        %ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info, residual] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_nn_ALM_nrep%dRowdim%dColdim%dRank%dTol%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        data(1, 1) = residual;
        data(2, 1) = cost;
        data(3, 1) = time;
    end

    if specifier.ind(2)
        %LQH
        fprintf('Starting LQH \n');
        timetic = tic();
        [xfinal, info, residual] = exactpenaltyViaSmoothinglqh(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_nn_LQH_nrep%dRowdim%dColdim%dRank%dTol%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres);
        struct2csv(info, filename);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        data(1, 2) = residual;
        data(2, 2) = cost;
        data(3, 2) = time;
    end
    
    if specifier.ind(3)
        %LSE
        fprintf('Starting LSE \n');
        timetic = tic();
        [xfinal, info, residual] = exactpenaltyViaSmoothinglse(problem, x0, options);
        time = toc(timetic);
        filename = sprintf('RC_nn_LSE_nrep%dRowdim%dColdim%dRank%dTol%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
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
        filename = sprintf('RC_nn_Riemannian_SQP_nrep%dRowdim%dColdim%dRank%dTol%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres);
        struct2csv(info, filename);        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        rankflag = check_rank(xfinal);
        if rankflag ~= 1
            residual = NaN;
        end
        data(1, 4) = residual;
        data(2, 4) = cost;
        data(3, 4) = time;
    end
    
    filename = sprintf('RC_nn_Info_nrep%dRowdim%dColdim%dRank%dTol%d.csv',setting.repeat,setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres);
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
end