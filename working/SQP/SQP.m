function [xfinal, costfinal, info, options] = SQP(problem0, x0, options)
% Sequential Quadratic Programming solver for smooth objective functions 
% on Riemannian manifolds.
%
% function [x, cost, info, options] = sqponmani(problem0)
% function [x, cost, info, options] = sqponmani(problem0, x0)
% function [x, cost, info, options] = sqponmani(problem0, x0, options)
% function [x, cost, info, options] = sqponmani(problem0, [], options)
%
% This is a Sequential Qudratic Programming solver for mixed constraints problems
% on Riemannian manifolds, which aims to minimize the cost function
% in the given problem structure with (in)equality constraints.
% It requires access to the gradient and the Hessian of the cost function
% and the constraints.
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below.
%
% The initial iterate is x0 if it is provided. Otherwise, a random point on
% the manifold is picked. To specify options whilst not specifying an
% initial iterate, give x0 as [] (the empty matrix).
%
% The two outputs 'x' and 'cost' are the last reached point on the manifold
% and its cost. 
% 
% The output 'info' is a struct-array which contains information about the
% iterations:
%   iter (integer)
%       The iteration number, or number of steps considered
%       (whether accepted or rejected). The initial guess is 0.
%	cost (double)
%       The corresponding cost value.
%	gradnorm (double)
%       The (Riemannian) norm of the gradient.
%	time (double)
%       The total elapsed time in seconds to reach the corresponding cost.
%	stepsize (double)
%       The size of the step from the previous to the new iterate.
%   accepted(bool)
%       The feasibility of the current point.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached at each iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tolgradnorm (1e-7)
%       The algorithm terminates if the norm of the gradient of the
%       Lagrangian (Notice here!)
%       drops below this. For well-scaled problems, a rule of thumb is that you can
%       expect to reduce the gradient norm by 7 orders of magnitude
%       (sqrt(eps)) compared to the gradient norm at a "typical" point (a
%       rough initial iterate for example). Further decrease is sometimes
%       possible, but inexact floating point arithmetic will eventually
%       limit the final accuracy. If tolgradnorm is set too low, the
%       algorithm may end up iterating forever (or at least until another
%       stopping criterion triggers).
%   maxiter (1000)
%       The algorithm terminates if maxiter iterations were executed.
%   maxtime (3600)
%       The algorithm terminates if maxtime seconds elapsed.
%   minstepsize (1e-6)
%       The minimum norm of the tangent vector that points from the current
%       point to the next point. If the norm is less than minstepsize, the 
%       program will terminate.
%   statsfun (none)
%       Function handle to a function that will be called after each
%       iteration to provide the opportunity to log additional statistics.
%       They will be returned in the info struct. See the generic Manopt
%       documentation about solvers for further information. statsfun is
%       called with the point x that was reached last, after the
%       accept/reject decision. See comment below.
%   stopfun (none)
%       Function handle to a function that will be called at each iteration
%       to provide the opportunity to specify additional stopping criteria.
%       See the generic Manopt documentation about solvers for further
%       information.
%   verbosity (2)
%       Integer number used to tune the amount of output the algorithm
%       generates during execution (mostly as text in the command window).
%       The higher, the more output. 0 means silent. 3 and above includes a
%       display of the options structure at the beginning of the execution.
%   debug (false)
%       Set to true to allow the algorithm to perform additional
%       computations for debugging purposes. If a debugging test fails, you
%       will be informed of it, usually via the command window. Be aware
%       that these additional computations appear in the algorithm timings
%       too, and may interfere with operations such as counting the number
%       of cost , etc. (the debug calls get storedb too).
%   storedepth (10)
%       Maximum number of different points x of the manifold for which a
%       store structure will be kept in memory in the storedb. If the
%       caching features of Manopt are not used, this is irrelevant. If
%       memory usage is an issue, you may try to lower this number.
%       Profiling may then help to investigate if a performance hit was
%       incurred as a result.
%
%
% Please cite the Manopt paper as well as the research paper:
% @InBook{Obara2020,
%   title     = {},
%   author    = {},
%   year      = {},
%   publisher = {},
%   editor    = {},
%   address   = {},
%   booktitle = {},
%   pages     = {},
%   doi       = {}
% }
%
%
% Original author: Mitsuaki Obara, January 20, 2020.
% Contributors: 
% Change log: 
%           January, 20, 2020: forked from Changshuo Liu's rlbfgs.m
%           deleted codes on localdefaults.tolgradnorm and 
%           localdefaults.strict_inc_func 

    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem0)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    if ~canGetGradient(problem0) && ~canGetApproxGradient(problem0)
        % Note: we do not give a warning if an approximate gradient is
        % explicitly given in the problem description, as in that case the user
        % seems to be aware of the issue.
        warning('manopt:getGradient:approx', ...
           ['No gradient provided. Using an FD approximation instead (slow).\n' ...
            'It may be necessary to increase options.tolgradnorm.\n' ...
            'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
        problem0.approxgrad = approxgradientFD(problem0);
   end
    
    % If the struct 'problem0' does not have a condet field (which is an
    % expected situation), add it to the problem0 here.
    if ~isfield(problem0, 'condet')
        problem0.condet = constraintsdetail(problem0);
    end
    
    % Set localdefaults, a struct to be combined with argument options for
    % declaring hyperparameters.
    localdefaults.maxiter = 1000;
    localdefaults.maxtime = 3600;
    localdefaults.tolqpnorm = 1e-6;
    localdefaults.tolstepnorm = 1e-6;
    localdefaults.tolgradnorm = 1e-7;
    % Hyperparameters for StoreDB
    localdefaults.storedepth = 3;
    % The way to treat the hessian matrix to be positive semidefinete.
    localdefaults.trimhessian = 'mineigval_manopt';
    % Initial parameters for the merit function and the Lagrangian
    localdefaults.tau = 0.8;  % TODO: should find an appropriate value as long as tau > 0
    localdefaults.rho = 1;  % TODO: should find an appropriate value as long as rho > 0
    localdefaults.beta = 0.5;  % TODO: should find an appropriate value as long as 1 > beta > 0
    localdefaults.gamma = 0.5; % TODO: should find an appropriate value as long as 1 > gamma > 0  
    localdefaults.mus = ones(problem0.condet.n_ineq_constraint_cost, 1);
    localdefaults.lambdas = ones(problem0.condet.n_eq_constraint_cost, 1);    
    localdefaults.ls_max_steps  = 30;
    localdefaults.regularhesseigval = 1e-3;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % Create a random starting point if no starting point is provided.
    if ~exist('x0', 'var')|| isempty(x0)
        xCur = problem0.M.rand(); 
    else
        xCur = x0;
    end
    
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    % create some initial variables which will be used in the following
    % loop.
    mus = options.mus; % initinal mus and lambdas for the Lagrangian
    lambdas = options.lambdas;
    rho = options.rho; % initinal rho for merit function
    stepsize = 1; % initial stepsize for linesearch
    % lsstats = []; % Line-search stastics for recording in info.
    
    % For the initial savestats, declare some variables
    iter = 0;
    timetic = tic();
    xCurCost = getCost(problem0, xCur, storedb, key);
    xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
    xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % stop flag, finally it should be true
    stop = false;
    totaltime = tic();
    
    % Main loop where we solve subproblems iteratively
    while true
        
        if options.verbosity >= 2
            fprintf('Iter: %d, Cost: %f, LagGradNorm: %f \n', iter, xCurCost, xCurLagGradNorm);
        end

        iter = iter + 1;
        
        timetic = tic();

        % Get current Hessian and gradient of the cost function.
        % Also, make a "qpinfo" structure stading for the subproblem
        % at the current point.
        costLag = @(X) costLagrangian(X, mus, lambdas); % though we don't use here.
        gradLag = @(X) gradLagrangian(X, mus, lambdas); % in the tangent space
        hessLag = @(X, d) hessLagrangian(X, d, mus, lambdas); % in the tangent space

        auxproblem.M = problem0.M;
        auxproblem.cost = costLag;
        auxproblem.grad = gradLag;
        auxproblem.hess = hessLag;

        % make H and basis
        [H,basis] = hessianmatrix(auxproblem, xCur);
        qpinfo.H = H;
        qpinfo.basis = basis;
        qpinfo.n = numel(basis);

        % make f
        f = zeros(qpinfo.n, 1);
        for fidx =1:qpinfo.n
            f(fidx) = auxproblem.M.inner(xCur, auxproblem.grad(xCur), basis{fidx});
        end
        qpinfo.f = f;

        % make inequality constraints
        if problem0.condet.has_ineq_cost
            row = problem0.condet.n_ineq_constraint_cost;
            col = numel(basis);
            A = zeros(row, col);
            b = zeros(row, 1);
            for ineqrow = 1:row
                costhandle = problem0.ineq_constraint_cost{ineqrow};
                b(ineqrow) = - costhandle(xCur);
                gradhandle = problem0.ineq_constraint_grad{ineqrow};
                constraint_egrad = gradhandle(xCur);
                constraint_grad = problem0.M.egrad2rgrad(xCur, constraint_egrad);
                for ineqcol = 1:col
                    base = basis{ineqcol};
                    A(ineqrow,ineqcol) = problem0.M.inner(xCur, constraint_grad, base);
                end
            end
        else
            A = [];
            b = [];
        end
        qpinfo.A = A;
        qpinfo.b = b;

        % make equality constraints
        if problem0.condet.has_eq_cost
            row = problem0.condet.n_eq_constraint_cost;
            col = numel(basis);
            Aeq = zeros(row, col);
            beq = zeros(row, 1);
            for eqrow = 1:row
                costhandle = problem0.eq_constraint_cost{eqrow};
                beq(eqrow) = - costhandle(xCur);
                gradhandle = problem0.eq_constraint_grad{eqrow};
                constraint_egrad = gradhandle(xCur);
                constraint_grad = problem0.M.egrad2rgrad(xCur, constraint_egrad);
                for eqcol = 1:col
                    base = basis{eqcol};
                    Aeq(eqrow,eqcol) = problem0.M.inner(xCur, constraint_grad, base);
                end
            end
        else
            Aeq = [];
            beq = [];
        end
        qpinfo.Aeq = Aeq;
        qpinfo.beq = beq;
        
        % Trim qpinfo.H (Hessian matrix) in some way, say, regularizeing
        % according to the minimum eigenvalue of Hessian obtained in some way
        % or replacing it with the identity matrix
        if strcmp(options.trimhessian, "eye")
            qpinfo.H = eye(qpinfo.n);
        elseif strcmp(options.trimhessian, 'mineigval_matlab')
            eigval = eig(qpinfo.H);
            qpinfo.mineigval = min(eigval);
            if qpinfo.mineigval < 0
                qpinfo.regularmineigval= max(options.regularhesseigval, abs(qpinfo.mineigval));
                qpinfo.H = qpinfo.H + qpinfo.regularmineigval * eye(qpinfo.n);
            end
        elseif strcmp(options.trimhessian, 'mineigval_manopt')
            % the difference between mineiegval_manopt and mineigval_matlab is
            % a solver for calculating the minimum eigen value. _matlab use the
            % eigenvalue decomposition function on matlab, whereas, _manopt use
            % the hessianextreme function which formulzize and solve
            % an optimization problem to get a minimum eigenvalue and 
            % the corresponding eigenvector.
            [~ ,qpinfo.mineigval] = hessianextreme(auxproblem, xCur);
            if qpinfo.mineigval < 0
                qpinfo.regularmineigval = max(options.regularhesseigval, abs(qpinfo.mineigval));
                qpinfo.H = qpinfo.H + qpinfo.regularmineigval * eye(qpinfo.n);
            end
        end

        % Compute the direction and Lagrange multipliers
        % by solving QP with quadprog, a matlab solver for QP
        [coeff, ~, ~, ~, Lagmultipliers] = quadprog(qpinfo.H, qpinfo.f,...
            qpinfo.A, qpinfo.b, qpinfo.Aeq, qpinfo.beq);
        
        deltaXast = 0;
        for i = 1:qpinfo.n
            deltaXast = deltaXast + coeff(i)* qpinfo.basis{i};
        end
        
        % Update rho, a penalty parameter, if needed.
        newacc = 0;
        for iterineq = 1 : problem0.condet.n_ineq_constraint_cost
            newacc = max(newacc, Lagmultipliers.ineqlin(iterineq));
        end
        
        for itereq = 1 : problem0.condet.n_eq_constraint_cost
            newacc = max(newacc, abs(Lagmultipliers.eqlin(itereq)));
        end
        
        if rho < newacc
           rho = newacc;
        end
        
        % make a struct and some variables for loneArmijoLineSearch
        meritproblem.M = problem0.M;
        meritproblem.cost = @(x) loneMeritFunction(x, rho);
        f0 = meritproblem.cost(xCur);
        
        % decide the value of df0 according to the way to
        % options.trimhessian
        if strcmp(options.trimhessian, "eye")
            df0 = meritproblem.M.inner(xCur, deltaXast, deltaXast);
        elseif strcmp(options.trimhessian, 'mineigval_matlab') || strcmp(options.trimhessian, 'mineigval_manopt')
            if qpinfo.mineigval < 0
                df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas) + qpinfo.regularmineigval * deltaXast, deltaXast);
            else 
                df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas), deltaXast);
            end
        else % standard setting
            df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas), deltaXast);
        end
        
        % Compute the stepsize with the L1-type merit function and the Armijo rule
        newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
        newf = meritproblem.cost(newx);
        gammadf0 = options.gamma * df0;
        r = 0;
        lsmaxiterbreak = false;
        % descriptCost(meritproblem, xCur, deltaXast); % DEBUG only
        while newf > ( f0 - gammadf0) && abs(newf - ( f0 - gammadf0)) > 10^(-4)
            if r > options.ls_max_steps
                lsmaxiterbreak = true;
                break;
            end
            r = r + 1;
            stepsize = stepsize * options.beta;
            gammadf0 = stepsize * gammadf0;
            newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
            newf = meritproblem.cost(newx);
        end
        
        % For savestats
        qpsolnorm = stepsize * problem0.M.norm(xCur, deltaXast);
        
        % Update variables to new iterate
        xCur = newx;
        mus = Lagmultipliers.ineqlin;
        lambdas =  Lagmultipliers.eqlin;
        
        % For savestats
        xCurCost = getCost(problem0, xCur, storedb, key);
        xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
        xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);        

        % savestats
        key = storedb.getNewKey();
        stats = savestats();
        info(iter+1) = stats;
        
        % refer to stop criteria        
        if iter >= options.maxiter
            fprintf('Max iter count reached\n');
            options.reason = "Max iter count reached";
            stop = true;
        elseif toc(totaltime) >= options.maxtime
            fprintf('Max time exceeded\n');
            options.reason = "Max time exceeded";
            stop = true;
        elseif xCurLagGradNorm <= options.tolgradnorm
            fprintf('Legrangian Gradient norm tolerance reached \n');
            options.reason = "Legrangian Gradient norm tolerance reached";
            stop = true;
        elseif qpsolnorm <= options.tolgradnorm
            fprintf('QP solution norm with stepsize tolerance reached\n');
            options.reason = "Legrangian Gradient norm tolerance reached";
            stop = true;
        end
        
        if stop
        options.totaltime = toc(totaltime);
        break
        end
        
    end
    
    xfinal = xCur;
    costfinal = problem0.cost(xfinal);
    
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = xCurCost;
        stats.gradnorm = xCurLagGradNorm;
        if iter == 0
            stats.time = toc(timetic);
            stats.qpsolnorm = NaN;
            stats.stepsize = NaN;
            stats.lsmaxiterbreak = NaN;
        else
            stats.time = toc(timetic);
            % stats.time = info(iter).time + toc(timetic);
            stats.qpsolnorm = qpsolnorm;
            stats.stepsize = stepsize;
            stats.lsmaxiterbreak = lsmaxiterbreak;
        end
        stats = applyStatsfun(problem0, xCur, storedb, key, options, stats);
    end

    function val = costLagrangian(x, mus, lambdas)
        val = getCost(problem0, x);
        if problem0.condet.has_ineq_cost
            for numineq = 1: problem0.condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + mus(numineq) * cost_numineq;
            end
        end

        if problem0.condet.has_eq_cost
            for numeq = 1: problem0.condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + lambdas(numeq) * cost_numeq;
            end
        end
    end

    function gradLag = gradLagrangian(x, mus, lambdas)
        gradLag = getGradient(problem0, x);
        if problem0.condet.has_ineq_cost
            for numineq = 1: problem0.condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                gradLag = problem0.M.lincomb(x, 1, gradLag, mus(numineq), constraint_grad);
            end
        end

        if problem0.condet.has_eq_cost
            for numeq = 1:problem0.condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                gradLag = problem0.M.lincomb(x, 1, gradLag, lambdas(numeq), constraint_grad);
            end
        end
    end

    function hessLag = hessLagrangian(x, dir, mus, lambdas)
        hessLag = getHessian(problem0, x, dir);
        if problem0.condet.has_ineq_cost
            for numineq = 1 : problem0.condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_egrad = gradhandle(x); % to be refactored
                hesshandle = problem0.ineq_constraint_hess{numineq};
                constraint_ehess = hesshandle(x, dir);
                constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                         constraint_ehess, dir);
                hessLag = problem0.M.lincomb(x, 1, hessLag,...
                    mus(numineq), constraint_hess);
            end
        end
        if problem0.condet.has_eq_cost
            for numeq = 1 : problem0.condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_egrad = gradhandle(x);
                hesshandle = problem0.eq_constraint_hess{numeq};
                constraint_ehess = hesshandle(x, dir);
                constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                        constraint_ehess, dir);
                hessLag = problem0.M.lincomb(x, 1, hessLag,...
                    lambdas(numeq), constraint_hess);
            end
        end
    end

    function val = loneMeritFunction(x, rho)
        val = getCost(problem0, x);
        if problem0.condet.has_ineq_cost
            for numineq = 1: problem0.condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + rho * max(0, cost_numineq);
            end
        end

        if problem0.condet.has_eq_cost
            for numeq = 1: problem0.condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + rho * abs(cost_numeq);
            end
        end
    end
end