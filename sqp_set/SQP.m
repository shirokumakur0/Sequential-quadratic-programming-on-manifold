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
%            LATER!
%
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached at each iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%                     LATER!
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

    % % Verify that the problem description is sufficient for the solver.
    % if ~canGetCost(problem0)
    %     warning('manopt:getCost', ...
    %         'No cost provided. The algorithm will likely abort.');
    % end
    % if ~canGetGradient(problem0) && ~canGetApproxGradient(problem0)
    %     % Note: we do not give a warning if an approximate gradient is
    %     % explicitly given in the problem description, as in that case the user
    %     % seems to be aware of the issue.
    %     warning('manopt:getGradient:approx', ...
    %        ['No gradient provided. Using an FD approximation instead (slow).\n' ...
    %         'It may be necessary to increase options.tolgradnorm.\n' ...
    %         'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
    %     problem0.approxgrad = approxgradientFD(problem0);
    % end
        
    condet = constraintsdetail(problem0);
    
    % Set localdefaults, a struct to be combined with argument options for
    % declaring hyperparameters.
    
    % For stopping criteria
    localdefaults.maxiter = 5000;
    localdefaults.maxtime = 3600;
    localdefaults.tolqpnorm = 1e-5;
    localdefaults.tolgradnorm = 1e-5;
    localdefaults.toliterdist = 1e-5;
    % For StoreDB
    localdefaults.storedepth = 3;
    % For modification for the hessian matrix to be positive semidefinete.
    localdefaults.modify_hessian = 'mineigval_matlab';
    localdefaults.mineigval_correction = 1e-8; % the param for mineigval_manopt or _matlab
    localdefaults.mineigval_threshold = 1e-3;
    % Initial parameters for the merit function and the Lagrangian
    localdefaults.tau = 0.8;  % TODO: should find an appropriate value as long as tau > 0
    localdefaults.rho = 1;  % TODO: should find an appropriate value as long as rho > 0
    localdefaults.beta = 0.5;  % TODO: should find an appropriate value as long as 1 > beta > 0
    localdefaults.gamma = 0.4; % TODO: should find an appropriate value as long as 1 > gamma > 0  
    localdefaults.mus = ones(condet.n_ineq_constraint_cost, 1);
    localdefaults.lambdas = ones(condet.n_eq_constraint_cost, 1);    
    % For linesearch
    localdefaults.ls_max_steps  = 50;
    localdefaults.ls_threshold = 1e-4;
    % For display
    localdefaults.verbosity = 1;
    localdefaults.qp_verbosity = 0;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % Set the quadprog verbosities
    if options.qp_verbosity == 0
        qpoptions = optimset('Display','off');
    else
        qpoptions = [];
    end
    
    % Create a random starting point if no starting point is provided.
    if ~exist('x0', 'var')|| isempty(x0)
        xCur = problem0.M.rand(); 
    else
        xCur = x0;
    end
    
    % Create a store database and get a key for the current x
    %storedb = StoreDB(options.storedepth);
    %key = storedb.getNewKey();
    
    % Create some initial variables which will be used in the following
    % loop.
    mus = options.mus; % Init. mus and lambdas for the Lagrangian
    lambdas = options.lambdas;
    rho = options.rho; % Init. rho for merit function
    
    % For the initial savestats, declare some variables
    iter = 0;
    xCurCost = getCost(problem0, xCur);
    xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
    xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
    timetic = tic();
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Stop flag, finally it should be true
    stop = false;
    totaltime = tic();
    
    % Main loop where we solve subproblems iteratively
    while true
        if options.verbosity >= 2
            fprintf('Iter: %d, Cost: %f, LagGradNorm: %f \n', iter, xCurCost, xCurLagGradNorm);
        elseif options.verbosity >= 1
            if mod(iter, 100) == 0 && iter ~= 0
                fprintf('Iter: %d, Cost: %f, LagGradNorm: %f \n', iter, xCurCost, xCurLagGradNorm);
            end
        end
        
        iter = iter + 1;
        timetic = tic();

        % Get current Hessian and gradient of the cost function.
        % Also, make a "qpinfo" structure stading for the subproblem
        % at the current point.
        
        costLag = @(X) costLagrangian(X, mus, lambdas); % which we'll use in hessextreme
        gradLag = @(X) gradLagrangian(X, mus, lambdas); % in the tangent space
        hessLag = @(X, d) hessLagrangian(X, d, mus, lambdas); % in the tangent space

        auxproblem.M = problem0.M;
        auxproblem.cost = costLag;
        auxproblem.grad = gradLag;
        auxproblem.hess = hessLag;

        % Make H, basis, and n
        [H,basis] = hessianmatrix(auxproblem, xCur);
        qpinfo.H = 0.5 * (H.'+H);
        qpinfo.basis = basis;
        qpinfo.n = numel(basis);

        % Make f
        f = zeros(qpinfo.n, 1);
        xCurGrad = getGradient(problem0, xCur);
        for fidx =1:qpinfo.n
            f(fidx) = problem0.M.inner(xCur, xCurGrad, basis{fidx});
        end
        qpinfo.f = f;

        % make inequality constraints
        if condet.has_ineq_cost
            row = condet.n_ineq_constraint_cost;
            col = qpinfo.n;
            A = zeros(row, col);
            b = zeros(row, 1);
            for ineqrow = 1:row
                ineqcosthandle = problem0.ineq_constraint_cost{ineqrow};
                b(ineqrow) = - ineqcosthandle(xCur);
                ineqgradhandle = problem0.ineq_constraint_grad{ineqrow};
                ineqconstraint_egrad = ineqgradhandle(xCur);
                ineqconstraint_grad = problem0.M.egrad2rgrad(xCur, ineqconstraint_egrad);
                for ineqcol = 1:col
                    base = qpinfo.basis{ineqcol};
                    A(ineqrow,ineqcol) = problem0.M.inner(xCur, ineqconstraint_grad, base);
                end
            end
        else
            A = [];
            b = [];
        end
        qpinfo.A = A;
        qpinfo.b = b;

        % make equality constraints
        if condet.has_eq_cost
            row = condet.n_eq_constraint_cost;
            col = qpinfo.n;
            Aeq = zeros(row, col);
            beq = zeros(row, 1);
            for eqrow = 1:row
                eqcosthandle = problem0.eq_constraint_cost{eqrow};
                beq(eqrow) = - eqcosthandle(xCur);
                eqgradhandle = problem0.eq_constraint_grad{eqrow};
                eqconstraint_egrad = eqgradhandle(xCur);
                eqconstraint_grad = problem0.M.egrad2rgrad(xCur, eqconstraint_egrad);
                for eqcol = 1:col
                    base = qpinfo.basis{eqcol};
                    Aeq(eqrow,eqcol) = problem0.M.inner(xCur, eqconstraint_grad, base);
                end
            end
        else
            Aeq = [];
            beq = [];
        end
        qpinfo.Aeq = Aeq;
        qpinfo.beq = beq;
        
        % Modify qpinfo.H (Hessian matrix) to be positive definite somehow.
        % The difference between mineiegval_manopt and mineigval_matlab is
        % a solver for calculating the minimum eigen value.
        if strcmp(options.modify_hessian, "eye")
            % The identity matrix as replacement to Hessian.
            qpinfo.H = eye(qpinfo.n);
        elseif strcmp(options.modify_hessian, 'mineigval_matlab') 
            % The eigenvalue decomposition function on matlab
            eigval = eig(qpinfo.H);
            qpinfo.mineigval = min(eigval);
            if qpinfo.mineigval < 0
                qpinfo.mineigval_diagcoeff= max(options.mineigval_threshold,...
                    abs(qpinfo.mineigval)) + options.mineigval_correction;
                qpinfo.H = qpinfo.H + qpinfo.mineigval_diagcoeff * eye(qpinfo.n);
            end
        elseif strcmp(options.modify_hessian, 'mineigval_manopt')
            % Rayleigh quotient minization to get get a minimum eigenvalue.
            [~ ,qpinfo.mineigval] = hessianextreme(auxproblem, xCur);
            if qpinfo.mineigval < 0
                qpinfo.mineigval_diagcoeff = max(options.mineigval_threshold,...
                    abs(qpinfo.mineigval)) + options.mineigval_correction;
                qpinfo.H = qpinfo.H + qpinfo.mineigval_diagcoeff * eye(qpinfo.n);
            end
        end
        qpinfo.H = 0.5 * (qpinfo.H.'+qpinfo.H);

        % Compute the direction and Lagrange multipliers
        % by solving QP with quadprog, a matlab solver for QP
        [coeff, ~, qpexitflag, ~, Lagmultipliers] = quadprog(qpinfo.H, qpinfo.f,...
            qpinfo.A, qpinfo.b, qpinfo.Aeq, qpinfo.beq, [], [], [], qpoptions);
        deltaXast = 0;
        for i = 1:qpinfo.n
            deltaXast = deltaXast + coeff(i)* qpinfo.basis{i};
        end
        
        % Update rho, a penalty parameter, if needed.
        newacc = rho;
        if condet.has_ineq_cost
            for iterineq = 1 : condet.n_ineq_constraint_cost
                newacc = max(newacc, Lagmultipliers.ineqlin(iterineq));
            end
        end
        if condet.has_eq_cost
            for itereq = 1 : condet.n_eq_constraint_cost
                newacc = max(newacc, abs(Lagmultipliers.eqlin(itereq)));
            end
        end
        if rho < newacc
           rho = newacc;
        end
        
        % Compute a problem and some variables for loneArmijoLineSearch
        meritproblem.M = problem0.M;
        meritproblem.cost = @(x) loneMeritFunction(x, rho);
        f0 = meritproblem.cost(xCur);
        
        % Compute df0 according to options.modify_hessian
        if strcmp(options.modify_hessian, "eye")
            df0 = meritproblem.M.inner(xCur, deltaXast, deltaXast);
        elseif strcmp(options.modify_hessian, 'mineigval_matlab') || strcmp(options.modify_hessian, 'mineigval_manopt')
            if qpinfo.mineigval < 0
                df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas) + qpinfo.mineigval_diagcoeff * deltaXast, deltaXast);
            else 
                df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas), deltaXast);
            end
        else % Standard setting
            df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
                    mus, lambdas), deltaXast);
        end
        
        % Compute the stepsize with the L1-type merit function and the Armijo rule
        stepsize = 1;
        newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
        newf = meritproblem.cost(newx);
        gammadf0 = df0 * options.gamma;
        r = 0; % back-tracking counter
        ls_max_steps_flag = false;
        
        % DEBUG only
        % descriptCost(meritproblem, xCur, deltaXast);
        
        while newf > ( f0 - gammadf0) && abs(newf - ( f0 - gammadf0)) > options.ls_threshold
            if r > options.ls_max_steps
                ls_max_steps_flag = true;
                break;
            end
            r = r + 1;
            stepsize = stepsize * options.beta;
            gammadf0 =  gammadf0 * options.beta;
            newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
            newf = meritproblem.cost(newx);
        end
        
        % For savestats
        % qpsolnorm = stepsize * problem0.M.norm(xCur, deltaXast);
        if contains(problem0.M.name(),'Stiefel')
            dist = norm(xCur - newx, 'fro');
        else
            dist = problem0.M.dist(xCur, newx);
        end
        % Update variables to new iterate
        xCur = newx;
        mus = Lagmultipliers.ineqlin;
        lambdas =  Lagmultipliers.eqlin;
        
        % For savestats (Cont'd)
        xCurCost = getCost(problem0, xCur);
        xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
        xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);        
        
        % savestats
        %key = storedb.getNewKey();
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
        %elseif qpsolnorm <= options.tolqpnorm
        %    fprintf('QP solution norm with stepsize tolerance reached\n');
        %    options.reason = "Legrangian Gradient norm tolerance reached";
        %    stop = true;
        elseif dist <= options.toliterdist
            fprintf('Distance of iteration tolerance reached\n');
            options.reason = 'Distance of iteration tolerance reached';
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
            % stats.qpsolnorm = NaN;
            stats.stepsize = NaN;
            stats.ls_max_steps_break = NaN;
            stats.dist =  NaN;
            stats.qpexitflag = NaN;
        else
            stats.time = toc(timetic);
            stats.time = info(iter).time + toc(timetic);
            % stats.qpsolnorm = qpsolnorm;
            stats.stepsize = stepsize;
            stats.ls_max_steps_break = ls_max_steps_flag;
            stats.dist = dist;
            stats.qpexitflag = qpexitflag;
        end
        stats.rho = rho;
        stats.violation_sum = violation_sum();
        [stats.maxviolation, stats.meanviolation] = const_evaluation(xCur);
        stats = applyStatsfun(problem0, xCur, [], [], options, stats);
    end

    function val = costLagrangian(x, mus, lambdas)
        val = getCost(problem0, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + mus(numineq) * cost_numineq;
            end
        end

        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + lambdas(numeq) * cost_numeq;
            end
        end
    end

    function gradLag = gradLagrangian(x, mus, lambdas)
        gradLag = getGradient(problem0, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                gradLag = problem0.M.lincomb(x, 1, gradLag, mus(numineq), constraint_grad);
            end
        end

        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                gradLag = problem0.M.lincomb(x, 1, gradLag, lambdas(numeq), constraint_grad);
            end
        end
    end

    function hessLag = hessLagrangian(x, dir, mus, lambdas)
        hessLag = getHessian(problem0, x, dir);
        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_egrad = gradhandle(x);
                hesshandle = problem0.ineq_constraint_hess{numineq};
                constraint_ehess = hesshandle(x, dir);
                constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                         constraint_ehess, dir);
                hessLag = problem0.M.lincomb(x, 1, hessLag,...
                    mus(numineq), constraint_hess);
            end
        end
        if condet.has_eq_cost
            for numeq = 1 : condet.n_eq_constraint_cost
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
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + rho * max(0, cost_numineq);
            end
        end

        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + rho * abs(cost_numeq);
            end
        end
    end

    function [maxviolation, meanviolation] = const_evaluation(xCur)
        maxviolation = 0;
        meanviolation = 0;
        
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                maxviolation = max(maxviolation, cost_at_x);
                meanviolation = meanviolation + max(0, cost_at_x);
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                maxviolation = max(maxviolation, cost_at_x);
                meanviolation = meanviolation + cost_at_x;
            end
        end
        if condet.has_ineq_cost || condet.has_eq_cost
            meanviolation = meanviolation / (condet.n_ineq_constraint_cost + condet.n_eq_constraint_cost);
        end
    end
    % For additiobal stats
    function val = violation_sum()
        xGrad = getGradient(problem0, xCur);
        val = problem0.M.norm(xCur, xGrad)^2;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                violation = max(0, cost_at_x);
                val = val + violation^2;
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                val = val + cost_at_x^2;
            end
        end
        val = sqrt(val);
        %stats.violation_sum = val;
    end
end