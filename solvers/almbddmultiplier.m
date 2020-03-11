function [xfinal, info, residual] = almbddmultiplier(problem0, x0, options)

    condet = constraintsdetail(problem0);
    
    %Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.lambdas = ones(condet.n_ineq_constraint_cost, 1);
    localdefaults.gammas = ones(condet.n_eq_constraint_cost, 1);
    localdefaults.bound = 20;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    localdefaults.outerverbosity = 1;  % verbosity for outer loops by MO
    localdefaults.tolKKTres = 1e-8; % A stopping criterion, added by MO
    localdefaults.minstepsize = 1e-8;  % added by MO
    %Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-6;
    
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
    lambdas = options.lambdas;
    gammas = options.gammas;
    rho = options.rho;
    oldacc = Inf;
    M = problem0.M;
    xCur = x0;
    xPrev = xCur;
    OuterIter = 0;
    
    % for savestats, by MO
    gradLagfun = @(X) gradLag(X, problem0, lambdas, gammas);
    xCurLagGrad = gradLagfun(xCur);
    xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
    xCurMaxLagMult = maxabsLagrangemultipliers(lambdas, gammas);
    xCurResidual = KKT_residual();
    
    stats = savestats(x0);
    info(1) = stats;
    info(min(10000, options.maxOuterIter+1)).iter = [];
    
    totaltime = tic();
    
    for OuterIter = 1 : options.maxOuterIter
        timetic = tic();
        % verbosity, odified by MO
        if options.outerverbosity >= 2
            fprintf('Iteration: %d    ', OuterIter);
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('Iteration: %d    ', OuterIter);
        end
        
        costfun = @(X) cost_alm(X, problem0, rho, lambdas, gammas);
        gradfun = @(X) grad_alm(X, problem0, rho, lambdas, gammas);
        problem.cost = costfun;
        problem.grad = gradfun;
        problem.M = M;
        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        inneroptions.minstepsize = options.minstepsize;
        
         
        [xCur, cost, innerinfo, Oldinneroptions] = rlbfgs(problem, xCur, inneroptions);
        
        % For KKT Residual, by MO
        gradLagfun = @(X) gradLag(X, problem0, lambdas, gammas);
        
        % updating judge for exit, added by MO
        updateflag_rho = false;
        updateflag_Lagmult = false;
        
        %Update Multipliers
        newacc = 0;
        for iterineq = 1 : condet.n_ineq_constraint_cost
            costhandler = problem0.ineq_constraint_cost{iterineq};
            cost_iter = costhandler(xCur);
            newacc = max(newacc, abs(max(-lambdas(iterineq)/rho, cost_iter)));
            
            % update judge, added by MO
            newlambda = min(options.bound, max(lambdas(iterineq) + rho * cost_iter, 0));
            if lambdas(iterineq) ~= newlambda
                lambdas(iterineq) = newlambda;
                updateflag_Lagmult = true;
            end
            % judge ended
        end
        
        for itereq = 1 : condet.n_eq_constraint_cost
            costhandler = problem0.eq_constraint_cost{itereq};
            cost_iter = costhandler(xCur);
            newacc = max(newacc, abs(cost_iter));
            
            % update judge, added by MO
            newgamma = min(options.bound, max(-options.bound, gammas(itereq) + rho * cost_iter));
            if gammas(itereq) ~= newgamma
                gammas(itereq) = newgamma;
                updateflag_Lagmult = true;
            end
            % judge ended
        end
        
        if OuterIter == 1 || newacc > options.tau * oldacc
            rho = rho/options.thetarho;
            updateflag_rho = true;
        end
        oldacc = newacc;
        oldtolgradnorm = tolgradnorm;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm);
        
        
        % For savestats, by MO
        xCurLagGrad = gradLagfun(xCur);
        xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
        xCurMaxLagMult = maxabsLagrangemultipliers(lambdas, gammas);
        xCurResidual = KKT_residual();
        
        % calculating the distance for savestats, by MO
        if contains(problem0.M.name(),'Stiefel') 
            dist = norm(xCur - xPrev, 'fro');
        elseif contains(problem0.M.name(),'rank')
            % Only assuming for 'fixedrankembeddedfactory'
            if ~exist('xCurmat', 'var')
                xCurmat = xCur.U * xCur.S * xCur.V';
            end    
            xPrevmat = xPrev.U * xPrev.S * xPrev.V';
            dist = norm(xCurmat - xPrevmat, 'fro');
            xCurmat = xPrevmat;
        else
            dist = problem0.M.dist(xCur, xPrev);
        end
        
        %Save stats
        stats = savestats(xCur);
        info(OuterIter+1) = stats;
        
        % verbosity modified by MO on March 2
        if options.outerverbosity >= 2
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        end
        
        % This is the stopping criterion based on violation_sum by MO on
        % March 4.
        if xCurResidual < options.tolKKTres && tolgradnorm <= options.endingtolgradnorm
            fprintf("KKT Residual tolerance reached\n")
            break;
        elseif dist == 0 && ~(updateflag_rho) && ~(updateflag_Lagmult) ...
                && (tolgradnorm == oldtolgradnorm)
            fprintf("Any parameter did not change\n")
            break; % because nothing changed, meaning that the alg. keeps producing the same point hereafter.
        end
        
        % The following part is as for a stopping criterion based on norm
        % by MO, March 2, 
        % if contains(problem0.M.name(),'rank')  % Only assuming for 'fixedrankembeddedfactory'.
        %     if ~exist('xPrevMat', 'var')
        %         xPrevMat = xPrev.U * xPrev.S * xPrev.V';
        %     end
        %     xCurMat = xCur.U * xCur.S * xCur.V';
        %     if norm(xPrevMat - xCurMat, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %         break;
        %     end
        %     xPrevMat = xCurMat;
        % elseif norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %     break;
        % end
        
        % The original one, remained here, just in case
        % if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %     break;
        % end
        
        if toc(totaltime) > options.maxtime
            break;
        end
        
        xPrev = xCur;
    end
    info = info(1: OuterIter+1);
    
    residual  = KKT_residual(); % added by MO

    xfinal = xCur;

    function stats = savestats(x)
        stats.iter = OuterIter;
        if stats.iter == 0
            stats.time = 0;
            stats.dist = NaN; % by MO
        else
            stats.time = info(OuterIter).time + toc(timetic);
            stats.dist = dist;
        end
        stats.rho = rho;
        [maxviolation, meanviolation, costCur] = evaluation(problem0, x, condet);
        stats.maxviolation = maxviolation;
        stats.meanviolation = meanviolation;
        stats.cost = costCur;
        % addding the information on Lagrange multipliers for sqp (by MO)
        stats.maxabsLagMult = xCurMaxLagMult;  % added by MO
        stats.KKT_residual = xCurResidual;  % added by MO
        stats.LagGradNorm = xCurLagGradNorm;
    end
    
    function val = cost_alm(x, problem, rho, lambdas, gammas)
        val = getCost(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + (rho/2) * (max(0, lambdas(numineq)/rho + cost_numineq)^2);
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + (rho/2) * (gammas(numeq)/rho + cost_numeq)^2;
            end
        end
    end

    function val = grad_alm(x, problem, rho, lambdas, gammas)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                if (cost_numineq + lambdas(numineq)/rho > 0)
                    gradhandle = problem.ineq_constraint_grad{numineq};
                    constraint_grad = gradhandle(x);
                    constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                    val = problem.M.lincomb(x, 1, val, cost_numineq * rho + lambdas(numineq), constraint_grad);
                end
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                gradhandle = problem.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                val = problem.M.lincomb(x, 1, val, cost_numeq * rho + gammas(numeq), constraint_grad);
            end
        end
    end

    % For additiobal stats (MO)
     function val = KKT_residual()
        xGrad = gradLagfun(xCur);
        manvio = manifoldViolation(xCur);
        val = (problem0.M.norm(xCur, xGrad))^2;
        val = val + manvio^2;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                violation = max(0, cost_at_x);
                val = val + (violation)^2;
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                val = val + (cost_at_x)^2;                
            end
        end
        val = sqrt(val);
        %stats.violation_sum = val;
     end

     % added by MO
    function val = maxabsLagrangemultipliers(lambdas, gammas)
       val = -1; % meaning no constraints
       if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                val = max(val, abs(lambdas(numineq)));
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                val = max(val, abs(gammas(numeq)) );
            end
        end
    end

    % added by MO
    function val = gradLag(x, problem0, lambdas, gammas)
        val = getGradient(problem0, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                val = problem0.M.lincomb(x, 1, val, lambdas(numineq), constraint_grad);
            end
        end

        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);
                val = problem0.M.lincomb(x, 1, val, gammas(numeq), constraint_grad);
            end
        end
    end

    % added by MO, for calculating KKT residual
    function manvio = manifoldViolation(xCur)
        % According to the type of manifold, calculate the violation from
        % constraints seen as the manifold.
        if contains(problem0.M.name(),'Sphere')         
            y = xCur(:);
            manvio = abs(y.'*y - 1);
        elseif contains(problem0.M.name(),'Oblique')
            [~,N] = size(xCur);
            colones = ones(N, 1);
            manvio = max(abs(diag(xCur.'*xCur)-colones));
        else % including fixed-rank manifolds
            manvio = 0;
        end
    end
end