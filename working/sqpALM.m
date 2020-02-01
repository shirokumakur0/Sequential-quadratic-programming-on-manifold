function [xfinal, cost, info] = sqpALM(problem0, x0, options)
% TODO: function [x, cost, info,options] = sqponmani(problem0, x0, options)
% This is the sqp algorithm with augumented Lagrangian method.

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
% It requires access to the gradient and the Hessian of the cost function.
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
%   accepted (boolean)
%       1 if the current step is accepted in the cautious update. 0 otherwise
%   And possibly additional information logged by options.statsfun.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached at each iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tolgradnorm (1e-6)
%       The algorithm terminates if the norm of the gradient drops below
%       this. For well-scaled problems, a rule of thumb is that you can
%       expect to reduce the gradient norm by 8 orders of magnitude
%       (sqrt(eps)) compared to the gradient norm at a "typical" point (a
%       rough initial iterate for example). Further decrease is sometimes
%       possible, but inexact floating point arithmetic will eventually
%       limit the final accuracy. If tolgradnorm is set too low, the
%       algorithm may end up iterating forever (or at least until another
%       stopping criterion triggers).
%   maxiter (1000)
%       The algorithm terminates if maxiter iterations were executed.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%   minstepsize (1e-10)
%     The minimum norm of the tangent vector that points from the current
%     point to the next point. If the norm is less than minstepsize, the 
%     program will terminate.
%   linesearch (@linesearch_hint)
%       Function handle to a line search function. The options structure is
%       passed to the line search too, so you can pass it parameters. See
%       each line search's documentation for info.
%       By default, the intial multiplier tried is alpha = 1. This can be
%       changed with options.linesearch: see help of linesearch_hint.
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
%   storedepth (30)
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
    
    % If the struct 'problem0' does not have a condet field (and it is a
    % expected situation), add the field to the problem0 here.
    if ~isfield(problem0, 'condet')
        problem0.condet = constraintsdetail(problem0);
    end
    
    % Set localdefaults, a struct to be combined with argument options for
    % declaring hyperparameters.
    localdefaults = setLocalDefaults(problem0);
    
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
    
    % Up to  here, the codes are borrowed from manopt and preceding works.
    % Now, we added the followings for SQP on manifolds.
    
    % Get the canonical basis corresponding with the manifold. We use this
    % for vectorizing gradients and make Hessianmatrix, both of which are
    % applied with Riemannian metrics.
    % NOTICE: The function only assume the case that the tangent spaces
    % of M have the same basis as that of M. For more detail, please check 
    % makeCanonicalBasis.m.
    basis = makeCanonicalBasis(problem0);
  
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    % __Initialization of variables__
    % create some variables which will be used in the following loop.
    M = problem0.M; % to make a subproblem on the manifold.
    mus = options.mus; % initinal mus and lambdas for the Lagrangian
    lambdas = options.lambdas;
    rho = options.rho; % initinal rho for merit function
    stepsize = 1; % initial stepsize for linesearch
    % lsstats = []; % Line-search stastics for recording in info.
    
    % For the initial savestats, declare some variables
    iter = 0;
    timetic =tic();
    [xCurCost, xCurGradient] = getCostGrad(problem0, xCur, storedb, key);
    xCurGradNorm = problem0.M.norm(xCur, xCurGradient);
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxouteriter+1)).iter = [];
       
    totaltime = tic();
    
    % Main loop where we solve subproblems iteratively
    for iter = 1:options.maxouteriter

        %         if options.verbosity >= 2
        %             fprintf(' iter                   cost val            grad. norm           alpha\n');
        %         end

        timetic = tic();

        % Get current Hessian and gradient of the cost function.
        % Also, make a "problem" structure stading for the subproblem
        % at the current point.
        fprintf('Iteration: %d     ', iter);
        costLag = @(X) costLagrangian(X, problem0, mus, lambdas); % value
        gradLag = @(X) gradLagrangian(X, problem0, mus, lambdas); % in the tangent space
        hessLag = @(X, d) hessLagrangian(X, d, problem0, mus, lambdas); % in the tangent space
        problem.cost = costLag;
        problem.grad = gradLag;
        problem.hess = hessLag;
        problem.M = M;

        % Make the grad and hess of the Lagrangian at the current point
        gradLagvec = gradMetricVectorize(xCur, gradLag(xCur), problem, basis);
        hessLagmat = hessMatLagrangian(xCur, problem, basis);

        % Tailor linearized constraints to the subproblem from the original
        % problem0.
        [ineqconst_gradmat, ineqconst_costvec, ...
         eqconst_gradmat, eqconst_costvec] = gradConstraintMatrix(xCur, problem0,...
                                                                  basis);

         % Compute the direction and Lagrange multipliers
         % by solving QP with quadprog, a matlab solver for QP
        [deltaXast, fval, ~, ~, Lagmultipliers] = quadprog(hessLagmat, gradLagvec,...
         ineqconst_gradmat, -ineqconst_costvec, eqconst_gradmat, -eqconst_costvec,...
         [],[],problem0.M.zerovec(xCur));

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

        % make the struct 'meritproblem', which consists of the L1 merit function
        % as meritproblem.cost, M as meritproble.M (for M.retr), and rho,tau,,
        f0 = loneMeritFunction(problem0, xCur, rho);
        df0 = fval - problem.M.inner(xCur, problem.grad(xCur),deltaXast);

        % Compute the stepsize with the L1-type merit function and the Armijo
        % rule
        % TODO: the return values rom loneMeritArmijoLineSearch should be 
        % [stepsize, newx, newkey, lsstats] for speeding up.
        [stepsize, newx] = loneMeritArmijoLineSearch(problem0,rho,...
                                                            xCur,deltaXast,f0,df0,options);
        
        % Update variables to new iterate
        xPrev = xCur;
        xCur = newx;
        mus = Lagmultipliers.ineqlin;
        lambdas = Lagmultipliers.eqlin;
        [xCurCost, xCurGradient] = getCostGrad(problem0, xCur, storedb, key);
        xCurGradNorm = problem0.M.norm(xCur, xCurGradient);
        key = storedb.getNewKey();
        % save stats
        stats = savestats();
        info(iter+1) = stats;
                                                        
        % refer to stop criteria        
        if toc(totaltime) >= options.maxtime
            fprintf('Max time exceeded');
            break
        elseif stepsize <= options.minstepsize
            fprintf('Min stepsize exceeded');
            break
        elseif problem0.M.norm(xPrev, deltaXast) <= options.tolgradnorm
            fprintf('Tol grad norm exceeded');
            break
        end
    end
    
    xfinal = xCur;
    cost = problem0.cost(xfinal);
    
    
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.xCur = xCur;
        stats.cost = xCurCost;
        stats.gradnorm = xCurGradNorm;
        if iter == 0
            stats.stepsize = NaN;
            stats.time = toc(timetic);
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
        end
        % stats.linesearch = lsstats;
        stats = applyStatsfun(problem0, xCur, storedb, key, options, stats);
    end
end