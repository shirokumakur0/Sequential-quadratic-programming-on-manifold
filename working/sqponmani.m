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
    localdefaults.linesearch = @linesearch_hint;
    
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
    
    timetic = tic();
    
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    % __________Initialization of variables______________
    % Number of iterations since the last restart
    k = 0;  
    % Total number of SQP iterations
    iter = 0; 
    
    
    % SHOULD BE DELETED!!
    % This cell stores step vectors which point from x_{t} to x_{t+1} for t
    % indexing the last iterations, capped at options.memory.
    % That is, it stores up to options.memory of the most recent step
    % vectors. However, the implementation below does not need step vectors 
    % in their respective tangent spaces at x_{t}'s. Rather, it requires
    % them transported to the current point's tangent space by vector
    % tranport. For details regarding the requirements on the the vector
    % tranport, see the reference paper by Huang et al.
    % In this implementation, those step vectors are iteratively 
    % transported to the current point's tangent space after every
    % iteration. Thus, at every iteration, vectors in sHistory are in the
    % current point's tangent space.
    sHistory = cell(1, options.memory);
    
    % This cell stores the differences for latest t's of the gradient at
    % x_{t+1} and the gradient at x_{t}, transported to x_{t+1}'s tangent
    % space. The memory is also capped at options.memory.
    yHistory = cell(1, options.memory);
    
    % rhoHistory{t} stores the reciprocal of the inner product between
    % sHistory{t} and yHistory{t}.
    rhoHistory = cell(1, options.memory);
    
    % Scaling of direction given by getDirection for acceptable step
    alpha = 1; 
    
    % Scaling of initial matrix, Barzilai-Borwein.
    scaleFactor = 1;
    
    % Norm of the step
    stepsize = 1;
    
    % Boolean for whether the step is accepted by Cautious update check
    accepted = 1;
    
    % Query the cost function and its gradient
    [xCurCost, xCurGradient] = getCostGrad(problem, xCur, storedb, key);
    
    xCurGradNorm = M.norm(xCur, xCurGradient);
    
    % Line-search statistics for recording in info.
    lsstats = [];
    
    % Flag to control restarting scheme to avoid infinite loops (see below)
    ultimatum = false;
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    if options.verbosity >= 2
        fprintf(' iter                   cost val            grad. norm           alpha\n');
    end
    
    % Main iteration
    while true

        % Display iteration information
        if options.verbosity >= 2
        fprintf('%5d    %+.16e        %.8e      %.4e\n', ...
                iter, xCurCost, xCurGradNorm, alpha);
        end
        
        % Start timing this iteration
        timetic = tic();
        
        % Run standard stopping criterion checks
        [stop, reason] = stoppingcriterion(problem, xCur, options, ...
                                           info, iter+1);
        
        % If none triggered, run specific stopping criterion check
        if ~stop 
            if stats.stepsize < options.minstepsize
                % To avoid infinite loop and to push the search further
                % in case BFGS approximation of Hessian is off towards
                % the end, we erase the memory by setting k = 0;
                % In this way, it starts off like a steepest descent.
                % If even steepest descent does not work, then it is 
                % hopeless and we will terminate.
                if ~ultimatum
                    if options.verbosity >= 2
                        fprintf(['stepsize is too small, restarting ' ...
                            'the bfgs procedure at the current point.\n']);
                    end
                    k = 0;
                    ultimatum = true;
                else
                    stop = true;
                    reason = sprintf(['Last stepsize smaller than '  ...
                        'minimum allowed; options.minstepsize = %g.'], ...
                        options.minstepsize);
                end
            else
                % We are not in trouble: lift the ultimatum if it was on.
                ultimatum = false;
            end
        end  
        
        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break;
        end

        % SHOULD BE MODIFIED!
        % Compute BFGS direction
        p = getDirection(M, xCur, xCurGradient, sHistory,...
                yHistory, rhoHistory, scaleFactor, min(k, options.memory));

        % Execute line-search
        [stepsize, xNext, newkey, lsstats] = ...
            linesearch_hint(problem, xCur, p, xCurCost, ...
                            M.inner(xCur, xCurGradient, p), ...
                            options, storedb, key);
        
        % Record the BFGS step-multiplier alpha which as effectively
        % selected. Toward convergence, we hope to see alpha = 1.
        alpha = stepsize/M.norm(xCur, p);
        step = M.lincomb(xCur, alpha, p);
        
        
        % Query cost and gradient at the candidate new point.
        [xNextCost, xNextGrad] = getCostGrad(problem, xNext, storedb, newkey);
        
        % Compute sk and yk
        sk = M.transp(xCur, xNext, step);
        yk = M.lincomb(xNext, 1, xNextGrad, ...
                             -1, M.transp(xCur, xNext, xCurGradient));

        % To make the Hessian operator more stable by normalizing the 
        % step taken.
        norm_sk = M.norm(xNext, sk);
        sk = M.lincomb(xNext, 1/norm_sk, sk);
        yk = M.lincomb(xNext, 1/norm_sk, yk);
        
        inner_sk_yk = M.inner(xNext, sk, yk);
        inner_sk_sk = M.norm(xNext, sk)^2;    % ensures nonnegativity
        
        
        % If the cautious step is accepted (which is the intended
        % behavior), we record sk, yk and rhok and need to do some
        % housekeeping. If the cautious step is rejected, these are not
        % recorded.
        cap = options.strict_inc_func(xCurGradNorm);
        if inner_sk_sk ~= 0 && (inner_sk_yk / inner_sk_sk) >= cap
            
            accepted = 1;
            
            rhok = 1/inner_sk_yk;
            
            scaleFactor = inner_sk_yk / M.norm(xNext, yk)^2;
            
            % Time to store the vectors sk, yk and the scalar rhok.
            % Remember: we need to transport all vectors to the most
            % current tangent space.
            
            % If we are out of memory
            if k >= options.memory
                
                % sk and yk are saved from 1 to the end with the most 
                % current recorded to the rightmost hand side of the cells
                % that are occupied. When memory is full, do a shift so
                % that the rightmost is earliest and replace it with the
                % most recent sk, yk.
                for  i = 2 : options.memory
                    sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                    yHistory{i} = M.transp(xCur, xNext, yHistory{i});
                end
                if options.memory > 1
                    sHistory = sHistory([2:end, 1]);
                    yHistory = yHistory([2:end, 1]);
                    rhoHistory = rhoHistory([2:end 1]);
                end
                if options.memory > 0
                    sHistory{options.memory} = sk;
                    yHistory{options.memory} = yk;
                    rhoHistory{options.memory} = rhok;
                end
                
            % If we are not out of memory
            else
                
                for  i = 1:k
                    sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                    yHistory{i} = M.transp(xCur, xNext, yHistory{i});
                end
                sHistory{k+1} = sk;
                yHistory{k+1} = yk;
                rhoHistory{k+1} = rhok;
                
            end
            
            k = k + 1;
            
        % The cautious step is rejected: we do not store sk, yk, rhok but
        % we still need to transport stored vectors to the new tangent
        % space.
        else
            
            accepted = 0;
            
            for  i = 1 : min(k, options.memory)
                sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                yHistory{i} = M.transp(xCur, xNext, yHistory{i});
            end
            
        end
        
        % Update variables to new iterate
        iter = iter + 1;
        xCur = xNext;
        key = newkey;
        xCurGradient = xNextGrad;
        xCurGradNorm = M.norm(xNext, xNextGrad);
        xCurCost = xNextCost;
        
        
        % Make sure we don't use too much memory for the store database
        % (this is independent from the BFGS memory.)
        storedb.purge();
        
        
        % Log statistics for freshly executed iteration
        stats = savestats();
        info(iter+1) = stats; 
        
    end

    
    % Housekeeping before we return
    info = info(1:iter+1);
    x = xCur;
    cost = xCurCost;

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', ...
                info(end).time);
    end

    
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = xCurCost;
        stats.gradnorm = xCurGradNorm;
        if iter == 0
            stats.stepsize = NaN;
            stats.time = toc(timetic);
            stats.accepted = NaN;
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.accepted = accepted;
        end
        stats.linesearch = lsstats;
        stats = applyStatsfun(problem, xCur, storedb, key, options, stats);
    end

end




% BFGS step, see Wen's paper for details. This functon takes in a tangent
% vector g, and applies an approximate inverse Hessian P to it to get Pg.
% Then, -Pg is returned.
%
% Theory requires the vector transport to be isometric and to satisfy the
% locking condition (see paper), but these properties do not seem to be
% crucial in practice. If your manifold provides M.isotransp, it may be
% good to replace M.transp with M.isotransp. There are built in M.isotransp
% for spherefactory and obliquefactory.
%
% This implementation operates in the tangent space of the most recent
% point since all vectors in sHistory and yHistory have been transported
% there.
function dir = getDirection(M, xCur, xCurGradient, sHistory, yHistory, ...
                            rhoHistory, scaleFactor, k)
    
    q = xCurGradient;
    
    inner_s_q = zeros(1, k);
    
    for i = k : -1 : 1
        inner_s_q(1, i) = rhoHistory{i} * M.inner(xCur, sHistory{i}, q);
        q = M.lincomb(xCur, 1, q, -inner_s_q(1, i), yHistory{i});
    end
    
    r = M.lincomb(xCur, scaleFactor, q);
    
    for i = 1 : k
         omega = rhoHistory{i} * M.inner(xCur, yHistory{i}, r);
         r = M.lincomb(xCur, 1, r, inner_s_q(1, i)-omega, sHistory{i});
    end
    
    dir = M.lincomb(xCur, -1, r);

end