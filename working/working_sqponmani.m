function [x, cost, info, options] = sqponmani(problem, x0, options)
% Sequential Quadratic Programming solver for smooth objective functions 
% on Riemannian manifolds.
%
% function [x, cost, info, options] = sqponmani(problem)
% function [x, cost, info, options] = sqponmani(problem, x0)
% function [x, cost, info, options] = sqponmani(problem, x0, options)
% function [x, cost, info, options] = sqponmani(problem, [], options)
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

% Original author: Mitsuaki Obara, January 20, 2020.
% Contributors: 
% Change log: 
%           January, 20, 2020: forked from Changshuo Liu's rlbfgs.m
%           deleted codes on localdefaults.tolgradnorm and 
%           localdefaults.strict_inc_func 


    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    if ~canGetGradient(problem) && ~canGetApproxGradient(problem)
        % Note: we do not give a warning if an approximate gradient is
        % explicitly given in the problem description, as in that case the user
        % seems to be aware of the issue.
        warning('manopt:getGradient:approx', ...
           ['No gradient provided. Using an FD approximation instead (slow).\n' ...
            'It may be necessary to increase options.tolgradnorm.\n' ...
            'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
        problem.approxgrad = approxgradientFD(problem);
    end
    
    % Local defaults for the program
    localdefaults.minstepsize = 1e-10;
    localdefaults.maxiter = 1000;
    localdefaults.tolgradnorm = 1e-6;
    localdefaults.ls_max_steps  = 25;
    localdefaults.storedepth = 30;
    localdefaults.tau = 0.8;  % TODO: should find an appropriate value as long as tau > 0
    localdefaults.rho = 1;  % TODO: should find an appropriate value as long as rho > 0
    localdefaults.beta = 0.5;  % TODO: should find an appropriate value as long as 1 > beta > 0
    localdefaults.gamma = 0.5; % TODO: should find an appropriate value as long as 1 > gamma > 0  

    % TODO: reconsider below. is this if-else part needed?
    if ~canGetLinesearch(problem)
        localdefaults.linesearch = @linesearch;
    else
        localdefaults.linesearch = @linesearch_hint;
    end
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    M = problem.M;
    
    % Create a random starting point if no starting point is provided.
    if ~exist('x0', 'var')|| isempty(x0)
        xCur = M.rand(); 
    else
        xCur = x0;
    end
        
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    % Rename
    M = problem.M;
    tau = options.tau;
    rho = options.rho;
    beta = options.beta;
    gamma = options.gamma;

    timetic = tic();
    
    % __Initialization of variables__
    % Number of iterations since the last restart
    k = 0;  
    % Total number of SQP iterations
    iter = 0;
    % Norm of the step
    stepsize = 1;
    % Line-search stastics for recording in info
    lsstats = [];
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    if options.verbosity >= 2
        fprintf(' iter                   cost val            grad. norm           alpha\n');
    end
    
    timetic = tic();

    % Get current Hessian and gradient of the cost function
    fprintf('Iteration: %d     ', iter);
    costfun = @(X) costLagrangian(X, problem0, mus, lambdas);
    gradfun = @(X) gradLagrangian(X, problem0, mus, lambdas);
    problem.cost = costfun;
    problem.grad = gradfun;
    problem.M = M;
    inneroptions.tolgradnorm = tolgradnorm;
    inneroptions.verbosity = 0;
    inneroptions.maxiter = options.maxInnerIter;
    inneroptions.minstepsize = options.minstepsize;
         
    [xCur, cost, innerinfo, Oldinneroptions] = quadprog(problem, xCur, inneroptions);
    
    % Tailor the quadratic model of the Lagrangian and linearized constraints 
    % to meet the subproblem, a quadratic programming; QP, at the current point
    
    % Compute the direction by solving QP with quadprog, a matlab solver for QP
    
    % Compute the stepsize with the L1-type merit function and the Armijo
    % rule
    
    % Update variables to new iterate
    
    % costLagrangian???verify????????????????
    
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = xCurCost;
        stats.gradnorm = xCurGradNorm;
        if iter == 0
            stats.stepsize = NaN;
            stats.time = toc(timetic);
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
        end
        stats.linesearch = lsstats;
        stats = applyStatsfun(problem, xCur, storedb, key, options, stats);
    end
end


function val = costLagrangian(x, problem, mus, lambdas)
    val = getCost(problem, x);
    if condet.has_ineq_cost
        for numineq = 1: condet.n_ineq_constraint_cost
            costhandle = problem.ineq_constraint_cost{numineq};
            cost_numineq = costhandle(x)
            val = val + mus(numineq) * cost_numineq;
        end
    end

    if condet.has_eq_cost
        for numeq = 1: condet.n_eq_constraint_cost
            costhandle = problem.eq_constraint_cost{numeq};
            cost_numeq = costhandle(x);
            val = val + lambdas(numeq) * cost_numeq;
        end
    end
end


function val = gradLagrangian(x, problem, mus, lambdas)
    val = getGradient(problem, x);
    if condet.has_ineq_cost
        for numineq = 1: condet.n_ineq_constraint_cost
            % costhandle = problem.ineq_constraint_cost{numineq};
            % cost_numineq = costhandle(x);
            gradhandle = problem.ineq_constraint_grad{numineq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
            val = problem.M.lincomb(x, 1, val, mus(numineq), constraint_grad);
        end
    end

    if condet.has_eq_cost
        for numeq = 1:condet.n_eq_constraint_cost
            % costhandle = problem.eq_constraint_cost{numeq};
            % cost_numeq = costhandle(x);
            gradhandle = problem.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
            val = problem.M.lincomb(x, 1, val, lambdas(numeq), constraint_grad);
        end
    end
end


% TODO: uncomplete, to use hessianmatrix, write the alg. for making normal
% basis.
function val = hessMatLagrangian(x, problem, mus, lambdas)
    
    val = getHessian(problem, x);
    if condet.has_ineq_cost
        for numineq = 1: condet.n_ineq_constraint_cost
            % costhandle = problem.ineq_constraint_cost{numineq};
            % cost_numineq = costhandle(x);
            gradhandle = problem.ineq_constraint_grad{numineq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
            val = problem.M.lincomb(x, 1, val, mus(numineq), constraint_grad);
        end
    end

    if condet.has_eq_cost
        for numeq = 1:condet.n_eq_constraint_cost
            % costhandle = problem.eq_constraint_cost{numeq};
            % cost_numeq = costhandle(x);
            gradhandle = problem.eq_constraint_grad{numeq};
            constraint_grad = gradhandle(x);
            constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
            val = problem.M.lincomb(x, 1, val, lambdas(numeq), constraint_grad);
        end
    end
end

