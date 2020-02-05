function [xfinal, cost, info, options] = sqp(problem0, x0, options)
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
%   maxtime (3600)
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
    problem0 = checkDifferentiability(problem0);
    
    % If the struct 'problem0' does not have a condet field (which is an
    % expected situation), add it to the problem0 here.
    if ~isfield(problem0, 'condet')
        problem0.condet = constraintsdetail(problem0);
    end
    
    % Set localdefaults, a struct to be combined with argument options for
    % declaring hyperparameters.
    localdefaults = setLocalDefaults(problem0);
    
    % Merge global and local defaults, then merge w/ user options, if any.
    options = trimOptions(options, localdefaults);
    
    % Create a random starting point if no starting point is provided.
    if ~exist('x0', 'var')|| isempty(x0)
        xCur = problem0.M.rand(); 
    else
        xCur = x0;
    end
    
    % Up to  here, the codes except condet are borrowed from manopt and 
    % preceding works. Now, we added the followings for SQP on manifolds.
    
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
    [xCurCost, xCurGradient] = getCostGrad(problem0, xCur, storedb, key);
    xCurGradNorm = problem0.M.norm(xCur, xCurGradient);
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Set totaltime for a stop criterion
    totaltime = tic();
    
    % Main loop where we solve subproblems iteratively
    for iter = 1:options.maxiter

        if options.verbosity >= 2
            fprintf('Iteration: %d, Cost: %f', iter, xCurCost);
        end

        timetic = tic();

        % Get current Hessian and gradient of the cost function.
        % Also, make a "qpinfo" structure stading for the subproblem
        % at the current point.
        qpinfo = makeQPInfo(problem0, xCur, mus, lambdas);
        
        % Trim qpinfo.H (Hessian matrix) in some way, say, regularizeing
        % according to the minimum eigenvalue of Hessian or replacing it
        % with the identity matrix
        qpinfo = trimHessianMatrix(qpinfo, options);
        
        
        
        % compare eigenvalue
        

        % Compute the direction and Lagrange multipliers
        % by solving QP with quadprog, a matlab solver for QP
        [coeff, fval, ~, ~, Lagmultipliers] = quadprog(qpinfo.H, qpinfo.f,...
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

        % make the struct 'meritproblem', which consists of the L1 merit function
        % as meritproblem.cost, M as meritproble.M (for M.retr), and rho,tau,,
        f0 = loneMeritFunction(problem0, xCur, rho);
        df0 = 2 * options.gamma * (fval - problem0.M.inner(xCur,...
            gradLagrangian(xCur,problem0, mus, lambdas),deltaXast));

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
        lambdas =  Lagmultipliers.eqlin;
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
            stats.deltaXastnorm = NaN;
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.deltaXastnorm = problem0.M.norm(xCur, deltaXast);
        end
        % stats.linesearch = lsstats;
        stats = applyStatsfun(problem0, xCur, storedb, key, options, stats);
    end
end