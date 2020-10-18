function [xfinal,info, residual] = exactpenaltyViaSmoothinglqh (problem0, x0, options)

    condet = constraintsdetail(problem0);
    
    %Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    localdefaults.startingepsilon = 1e-1;
    localdefaults.endingepsilon = 1e-6;
    localdefaults.outerverbosity = 1;  % verbosity for outer loops by MO
    localdefaults.tolKKTres = 1e-8; % a stopping criterion, added by MO
    localdefaults.minstepsize = 1e-8;  % added by MO
    %Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-6;
    % For fixed-rank manifolds, rank check at KKT residual
    localdefaults.rankviopena = 1e+8;
    
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);

    % added for the KKT residual at fixed-rank manifolds 3.31
    if contains(problem0.M.name(),'rank')
        if isfield(options, 'rank')
            rankval = options.rank;
        else
            tmpx = problem0.M.rand();
            rankval = rank(tmpx.S);
        end
    end    
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    theta_epsilon = nthroot(options.endingepsilon/options.startingepsilon, options.numOuterItertgn);
    
    M = problem0.M;
    xCur = x0;
    xPrev = xCur;
    epsilon = options.startingepsilon;
    rho = options.rho;
    
    % for savestats, by MO
    xCurMaxLagMult = maxabsLagrangemultipliers(xCur, problem0, epsilon);
    gradfun = @(X) grad_exactpenalty(X, problem0, rho);
    xCurResidual = KKT_residual();
    
    
    OuterIter = 0;
    stats = savestats(x0);
    info(1) = stats;
    info(min(10000, options.maxOuterIter+1)).iter = [];

    totaltime = tic();
    
    for OuterIter = 1 : options.maxOuterIter
        timetic = tic();
        
        % verbosity modified by MO
        if options.outerverbosity >= 2
            fprintf('Iteration: %d    ', OuterIter);
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('Iteration: %d    ', OuterIter);
        end

        costfun = @(X) cost_exactpenalty(X, problem0, rho);
        gradfun = @(X) grad_exactpenalty(X, problem0, rho);
        problem.cost = costfun;
        problem.grad = gradfun;
        problem.M = M;
        
        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        inneroptions.minstepsize = options.minstepsize;
        
        [xCur, cost, innerInfo, Oldinneroptions] = rlbfgs(problem, xCur, inneroptions);
        
        % for savestats, by MO
        xCurMaxLagMult = maxabsLagrangemultipliers(xCur, problem0, epsilon);
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
        
        % According to the convergence theory. Liu and Boumal, 2019
        %if stats.maxviolation > epsilon
        %    rho = rho/options.thetarho;
        %end
        
        % updating judge for exit, added by MO
        oldeps = epsilon;
        oldtolgradnorm = tolgradnorm;
        
        epsilon  = max(options.endingepsilon, theta_epsilon * epsilon);
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm);
        
        % verbosity modified by MO on March 2
        if options.outerverbosity >= 2
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        end
        
        % This is the stopping criterion based on violation_sum by MO on
        % March 4.
        if xCurResidual < options.tolKKTres && tolgradnorm <= options.endingtolgradnorm
            fprintf("KKT residual tolerance reached\n")
            break;
        elseif dist == 0 && (epsilon == oldeps) && (tolgradnorm == oldtolgradnorm)
            fprintf("Any parameter did not change\n")
        %if dist == 0 && (epsilon == oldeps) && (tolgradnorm == oldtolgradnorm)
        %    fprintf("Any parameter did not change\n")
            break; % because nothing changed, meaning that the alg. keeps producing the same point hereafter.
        end
        % The following part (norm judge and breaking) is modifiied by MO,
        % March 2
        % for fixedrankembeddedfactory
        %if contains(problem0.M.name(),'rank')  % Only assuming for 'fixedrankembeddedfactory'.
        %    if ~exist('xPrevMat', 'var')
        %        xPrevMat = xPrev.U * xPrev.S * xPrev.V';
        %    end
        %    xCurMat = xCur.U * xCur.S * xCur.V';
        %    if norm(xPrevMat - xCurMat, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %        break;
        %    end
        %    xPrevMat = xCurMat;
        %elseif norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %    break;
        %end
        %
        % The original stopping criterion, remained here just in case
        % if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %     break;
        % end
        % modification is up to here.
        
        xPrev = xCur;
        
        if toc(totaltime) > options.maxtime
            break;
        end
        
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
            stats.dist = dist; % by MO
        end
        [maxviolation, meanviolation, costCur] = evaluation(problem0, x, condet);
        stats.maxviolation = maxviolation;
        stats.meanviolation = meanviolation;
        stats.cost = costCur;
        % addding the information on Lagrange multipliers for sqp (by MO)
        stats.maxabsLagMult = xCurMaxLagMult;  % added by MO
        stats.KKT_residual = xCurResidual;  % added by MO
    end
    

    function val = cost_exactpenalty(x, problem, rho)
        val = getCost(problem, x);
        % Adding ineq constraint cost
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(x);
                if cost_at_x <= 0
                    cost_at_x = 0;
                elseif cost_at_x > epsilon
                        cost_at_x = cost_at_x - epsilon/2;
                else
                    cost_at_x = cost_at_x^2/(2*epsilon);
                end
                val = val + rho*cost_at_x;
            end
        end
        %Eq constratint cost
        if condet.has_eq_cost
            for numeq = 1 : condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_at_x = costhandle(x);
                val = val + rho * sqrt(cost_at_x^2 + epsilon^2);
            end
        end
    end

    function val = grad_exactpenalty(x, problem, u)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                if cost_at_x >= 0
                    gradhandle = problem.ineq_constraint_grad{numineq};
                    constraint_grad = gradhandle(x);
                    constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                    if cost_at_x >= epsilon
                        coef = 1;
                    else
                        coef = cost_at_x/epsilon;
                    end
                    val = problem.M.lincomb(x, 1, val, u * coef, constraint_grad);
                end
            end
        end
        if condet.has_eq_cost
            for numineq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                gradhandle = problem.eq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                coef = cost_at_x/sqrt(cost_at_x^2+epsilon^2);
                val = problem.M.lincomb(x, 1, val, u * coef, constraint_grad);
            end 
        end
    end

    % Added by MO
    function val = KKT_residual()
        %grad = gradLag(xCur, problem0, epsilon);
        grad = gradfun(xCur);
        val = problem0.M.norm(xCur, grad)^2;
        manpowvio = manifoldPowerViolation(xCur);
        compowvio = complementaryPowerViolation(xCur, problem0, epsilon);
        val = val + manpowvio;
        val = val + compowvio;
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
        
        % check if the current point satisfies the rank const. alternative
        % role as manvio at fixed-rank manifolds
        if contains(problem0.M.name(),'rank')
            xCurMatRank = xCur.U * xCur.S * xCur.V';
            xCurrank = rank(xCurMatRank);
            if xCurrank ~= rankval
                val = options.rankviopena; % as if Inf;
            end
        end
        
    end
    
    function compowvio = complementaryPowerViolation(x, problem, u)
        compowvio = 0;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                lambda = 1 / ( 1 + exp(- cost_at_x / u) );
                
                if cost_at_x <= 0
                    lambda = 0;
                elseif cost_at_x <= u
                    lambda = cost_at_x / u;
                else
                    lambda = 1;
                end
                violation = rho * lambda * cost_at_x;
                compowvio = compowvio + violation^2;
            end
        end
    end

    % Added by MO
    function val = gradLag(x, problem, u)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                gradhandle = problem.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                if cost_at_x <= 0
                    lambda = 0;
                elseif cost_at_x <= u
                    lambda = cost_at_x / u;
                    val = problem.M.lincomb(x, 1, val, rho * lambda, constraint_grad);
                else
                    lambda = 1;
                    val = problem.M.lincomb(x, 1, val, rho * lambda, constraint_grad);
                end
            end
        end
        if condet.has_eq_cost
            for numineq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                gradhandle = problem.eq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                gamma = cost_at_x / sqrt(cost_at_x^2 + u^2);
                val = problem.M.lincomb(x, 1, val, rho * gamma, constraint_grad);
            end 
        end
    end

    % Added by MO
    function val = maxabsLagrangemultipliers(x, problem, u)
        val = -1; % meaning no constraints
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                if cost_at_x <= 0
                    lambda = 0;
                elseif cost_at_x <= u
                    lambda = cost_at_x / u;
                else
                    lambda = 1;
                end
                val = max(val, abs(lambda));
            end
        end
        if condet.has_eq_cost
            for numineq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                gamma = cost_at_x / sqrt(cost_at_x^2 + u^2);
                val = max(val, abs(gamma));
            end 
        end
    end

    % added by MO, for calculating KKT residual
    function manvio = manifoldPowerViolation(xCur)
        % According to the type of manifold, calculate the violation from
        % constraints seen as the manifold.
        manvio = 0;
        if contains(problem0.M.name(),'Sphere')         
            y = xCur(:);
            manvio = abs(y.'*y - 1)^2;
        elseif contains(problem0.M.name(),'Oblique')
            [~,N] = size(xCur);
            %colones = ones(N, 1);
            for i = 1:N
                manvio = manvio + abs(xCur(:,i).' * xCur(:,i) - 1)^2;
            end
            %manvio = max(abs(diag(xCur.'*xCur)-colones));
        else % including fixed-rank manifolds
            manvio = 0;
        end
    end
end