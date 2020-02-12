function varargout = RBFGS(fns,params)
% Solver of Riemannian BFGS method
%
% INPUT:
% fns : a struct that contains required function handles
% required:
%     fns.f(x) : return objective function value at x.
%     fns.gf(x) : return the gradient of objection function at x.
% not required if params.manifold, params.retraction and params.vector_transport is provided.
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
%     fns.Hv(x, H, v) : return tangent vector which is given by that operator H affect on tangent vector v on T_x M.
%     fns.STV(a, v) : return a tangent vector which is a times v.
%     fns.TVaddTV(v1, v2) : return a tangent vector which is equal to v1 adds v2.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
%     fns.beta(x, eta) : return \|eta\| / \|TR_eta (eta)\|, where TR is the differentiated retraction.
%     fns.Tranv(x1, d, x2, v) : return a tangent vector on T_x2 M, which is given by vector transport that transport v \in T_x M to T_{R_x(d)} M. Here x2 = R_x(d).
%     fns.TranH(x1, d, x2, H1) : return a operator, H2, that affects on T_x2 M. H2 satisfies that for any v \in T_x2(M), H2(v) = Tran (H1 (Tran^{-1}(v)))
%                             where Tran is the vector transport fns.Tranv(x1, d, x2, v).
%     fns.rank1operator(x, v1, v2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = g_x(v2, v) * v1.
%     fns.operadd(x, H1, H2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = H1(v) + H2(v).
%
% params : a struct that contains parameters that are used in the solver.
% required:
%     params.x0 : initial approximation of minimizer.
%     params.H0 : initial approximation of inverse Hessian.
% not required:
%     params.error [1e-5]       : tolerance of stopping criterion.
%     params.StopCriterion [3]  : stopping criterion, 0 means call fns.IsStopped to check
%                                 1 means stop when relative error of objective function is less than tolerance,
%                                 2 means stop when norm of gradient is less than tolerance,
%                                 3 means stop when norm of gradient over initial norm of gradient is less than tolerance,
%                                 4 means stopping criterion for partly smooth function.
%     params.restart [1]        : 1 means restart if curvature condition is not satisfied, 2 means return directly.
%                                 (This should not happen theoretically. In practice, only numerical error make it occur occasionally.)
%     params.max_t [300]        : the maximum number of iterations.
%     params.debug [1]          : '0' means almost silence, '1' means some output. '2' means a lot of output, '3' means more than need to know.
%     params.alpha [1e-4]       : the coefficient of Wolfe first condition. It is between (0, 0.5].
%     params.beta [0.999]       : the coefficient of Wolfe second condition. It is between [0.5, 1).
%     params.minstepsize [1e-10] : the min step size of line search algorithm. It is between [0, 0.1].
%     params.maxstepsize [200]  : the max step size of line search algorithm. It is between [2, Inf].
%     params.err_x [1e-5]       : tolerance for checking whether iterates are close to convergence. It is between [0, 1].
%     params.num_Grads [10]     : number of gradients for checking whether iterates are close to Clark stationary point. It is between [0, inf].
% The following give some examples of manifold. For each manifold, some retractions and vector transports are provided.
%     params.manifold : the manifold which objective function is on.
%         '1', sphere : S^{n - 1}
%         '2', Stiefel manifold : St(p, n)
%         '3', orthogonal group : O(n)
%         '4', grassmann manifold : Gr(p, n)
%
% OUTPUT:
% varargout : contains output parameters
%     varargout{1} : all iterates
%     varargout{2} : function values of all iterates
%     varargout{3} : norm of gradient of all iterates
%     varargout{4} : total computational time
%     varargout{5} : number of function evaluations
%     varargout{6} : number of gradient evaluations
%     varargout{7} : number of Hessian actions
%     varargout{8} : accumulated computational time for iterations
%     varargout{9} : number of vector transport actions
%     varargout{10} : number of retraction actions
%
% By Wen Huang
%fprintf('wocacacacacca\n')
fprintf('RBFGS\n')
if nargin < 2,
    error('Invalid arguments: the number of arguments is invalid.');
end
% set default parameters and check the legal range of parameters
params = set_default_para(params, 'StopCriterion', 3, 'int', 0, 3);
params = set_default_para(params, 'error', 1e-6, 'float', 0, 1);
params = set_default_para(params, 'err_x', 1e-5, 'float', 0, 1);
params = set_default_para(params, 'num_Grads', 10, 'int', 0, inf);
params = set_default_para(params, 'max_t', 300, 'int', 1, inf);
params = set_default_para(params, 'debug', 1, 'int', 0, 3);
params = set_default_para(params, 'alpha', 1e-4, 'float', 0, 0.5);
params = set_default_para(params, 'beta', 0.999, 'float', 0.5, 1);
params = set_default_para(params, 'minstepsize', 1e-10, 'float', 0, 0.1);
params = set_default_para(params, 'maxstepsize', 200, 'int', 2, inf);
params = set_default_para(params, 'isconvex', 0, 'int', 0, 1);
params = set_default_para(params, 'keepnumX', 2100000000, 'int', 1, inf);
params = set_default_para(params, 'isplot', 0, 'int', 0, 1);
params = set_default_para(params, 'initstepsize', 1, 'float', 0, inf);
params = set_default_para(params, 'linesearch', 1, 'int', 1, 2);
params = set_default_para(params, 'ratio', 0.25, 'float', 0, 1);
params = set_default_para(params, 'nu', 1e-6, 'float', 0, inf);
params = set_default_para(params, 'mu', 1, 'float', 0, inf);
if(isfield(params, 'manifold') && isfield(params, 'retraction') && isfield(params, 'vector_transport'))
    fs = choose_manifold(params.manifold, params.retraction, params.vector_transport);
    fns.Hv = fs.Hv;
    fns.STV = fs.STV;
    fns.TVaddTV = fs.TVaddTV;
    fns.R = fs.R;
    fns.beta = fs.beta;
    fns.inpro = fs.inpro;
    fns.Tranv = fs.Tranv;
    fns.TranH = fs.TranH;
    fns.rank1operator = fs.rank1operator;
    fns.operadd = fs.operadd;
end
% Check whether function handles are legal
check_function_handle(fns, 'f');
check_function_handle(fns, 'gf');
check_function_handle(fns, 'Hv');
check_function_handle(fns, 'STV');
check_function_handle(fns, 'TVaddTV');
check_function_handle(fns, 'R');
check_function_handle(fns, 'beta');
check_function_handle(fns, 'inpro');
check_function_handle(fns, 'Tranv');
check_function_handle(fns, 'TranH');
check_function_handle(fns, 'rank1operator');
check_function_handle(fns, 'operadd');
if(params.StopCriterion == 0)
    check_function_handle(fns, 'IsStopped');
end
% Check whether required parameters are legal
if(~isfield(params, 'x0') || ~isfield(params, 'H0'))
    error('Invalid arguments: missing initial x0 or H0');
end
% Initialization
err = inf;
err_v = inf;
err_x = inf;
err_g = inf;
times = 0;
x1 = params.x0;
H = params.H0;
[f1, x1] = fns.f(x1);
[gradf1, x1] = fns.gf(x1);

X{1} = x1;
F(1) = f1;
G(1) = sqrt(fns.inpro(x1, gradf1, gradf1));

nf = 1;
ng = 1;
nH = 0;
nV = 0;
nR = 0;
ind = 1;
Grads = [];
timeinfo = '';
fprintf('times : %d, f(x) : %e\n', times, f1);
tic
T(1) = 0;
while(err > params.error && times < params.max_t)
    if(params.debug > 1)
        getetatimestart = toc;
    end
    eta = fns.Hv(x1, H, gradf1);
    eta = fns.STV(x1, -1, eta);
    if(params.debug > 1)
        timeinfo = sprintf('time, get eta:%.2e,', toc - getetatimestart);
        linesearchtimestart = toc;
    end
    initslope = fns.inpro(x1, eta, gradf1);
    if(times == 0)
        initstepsize = params.initstepsize;
    else
        initstepsize = min([1, 1.01 * (2 * (f1 - fpre) / initslope)]);
    end
    if(params.StopCriterion == 4)
        [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, timeinfo] = linesearch_PartlySmooth(fns, eta, x1, f1, gradf1, initslope, initstepsize, params, timeinfo);
    else
        if(params.linesearch == 1)
            if(err_g < params.accuracy) 
                t = 1;
                x2 = fns.R(x1,eta);
                lnR = lnR + 1;
                [f2, x2] = fns.f(x2);
                lnf = lnf + 1;
                [gradf2, x2] = fns.gf(x2);
                lng = lng + 1;
                ls = 1;
            else
                [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, timeinfo] = linesearch_Armijo(fns, eta, x1, f1, gradf1, initslope, initstepsize, params, timeinfo);
                
            end
        else
            if(err_g < params.accuracy) 
                t = 1;
                x2 = fns.R(x1,eta);
                lnR = lnR + 1;
                [f2, x2] = fns.f(x2);
                lnf = lnf + 1;
                [gradf2, x2] = fns.gf(x2);
                lng = lng + 1;
                ls = 1;
            else
                [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, timeinfo] = linesearch_Wolfe(fns, eta, x1, f1, gradf1, initslope, initstepsize, params, timeinfo);
            end
            
        end
    end
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('linesearch total:%.2e,', toc - linesearchtimestart)];
        getsytimestart = toc;
    end
    if(ls < 0)
        fprintf('warning: line search fails!\n');
        break;
    end
    nV = nV + lnV;
    [s, eta, x1, x2] = fns.Tranv(x1, eta, x2, eta);
    [Tgf, eta, x1, x2] = fns.Tranv(x1, eta, x2, gradf1);
    if(params.debug > 2)
        cprintf([1,0,1], sprintf('\t\tCheck isometry of vector transport: |gf| / |T gf|: %e/%e \n', ...
            sqrt(fns.inpro(x1, gradf1, gradf1)), sqrt(fns.inpro(x2, Tgf, Tgf))));
    end
    y = fns.TVaddTV(x2, fns.STV(x2, 1 / fns.beta(x1, eta, x2), gradf2), fns.STV(x2, -1, Tgf));
    nV = nV + 2;
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('get s y:%.2e,', toc - getsytimestart)];
        transportHtimestart = toc;
    end
    inp_s_y = fns.inpro(x2, s, y);
    if(times == 0 && params.isconvex)
        H = (inp_s_y / fns.inpro(x2, y, y)) * H;
    end
    if(params.debug > 2)
        cprintf([1,0,1], sprintf('\t\tCheck isometry of inverse Hessian approximation: g(gradf1, H eta) / g(H gradf1, eta): %e/%e \n', ...
            fns.inpro(x1, gradf1, fns.Hv(x1, H, eta)), fns.inpro(x1, fns.Hv(x1, H, gradf1),eta)));
    end
    % update operator H
    tilde_H = fns.TranH(x1, eta, x2, H);
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('transport H:%.2e,', toc - transportHtimestart)];
        updateHtimestart = toc;
    end
    
    inp_s_s = fns.inpro(x2, s, s);
    ngf = sqrt(fns.inpro(x1, gradf1, gradf1));
    if(inp_s_y / inp_s_s >= params.nu * ngf^(params.mu))
        %         Hy = fns.Hv(x2, tilde_H, y);
        %         H = fns.operadd(x2, tilde_H, fns.rank1operator(x2, fns.STV(x2, -1 / inp_s_y, s), Hy));
        %         H = fns.operadd(x2, H, fns.rank1operator(x2, Hy, fns.STV(x2, - 1 / inp_s_y, s)));
        %         H = fns.operadd(x2, H, fns.rank1operator(x2, fns.STV(x2, 1 / inp_s_y, s), fns.STV(x2, (fns.inpro(x2, y, Hy) / inp_s_y + 1), s)));
        Hy = fns.Hv(x2, tilde_H, y);
        H = fns.operadd(x2, tilde_H, fns.rank1operator(x2, fns.STV(x2, -1 / inp_s_y, s), Hy));
        Hy = fns.Hv(x2, H, y);
        H = fns.operadd(x2, H, fns.rank1operator(x2, Hy, fns.STV(x2, - 1 / inp_s_y, s)));
        H = fns.operadd(x2, H, fns.rank1operator(x2, fns.STV(x2, 1 / inp_s_y, s), s));
    else
        H = tilde_H;
        if(params.debug > 0)
            fprintf('without update H\n');
        end
    end
    
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('update H:%.2e,', toc - updateHtimestart)];
    end
    % check stoping cretirion
    err_r = abs((f2 - f1) / f2);
    if(params.StopCriterion == 4)
        if(params.debug > 1)
            checkstoppartlysmooth = toc;
        end
        if(iscell(x1))
            err_x = norm(x1{1} - x2{1});
        else
            err_x = sqrt(fns.inpro(x1, eta, eta));% norm(x1 - x2); %%---Whether err_x is OK is still a problem.
        end                                       %% The ideal form should be dist(x1, x2) which may be expensive.
        for i = 1 : size(Grads, 2)
            if(i ~= ind)
                Gradsi.TV = reshape(Grads(:, i), size(gradf1.TV));
                [tempv, eta, x1, x2] = fns.Tranv(x1, eta, x2, Gradsi);
                Grads(:, i) = reshape(tempv.TV, [], 1);
                nV = nV + 1;
            end
        end
        
        if(err_x < params.err_x)
            Grads(:, ind) = reshape(gradf2.TV, [], 1);
            ind = ind + 1;
            if(ind > params.num_Grads)
                ind = ind - params.num_Grads;
            end
            [w, err_vectorv] = qpsubprob(Grads, 1);
            err_vector.TV = reshape(err_vectorv, size(gradf1.TV));
            %                 gradf2 = - err_vector;
            err_v = sqrt(fns.inpro(x2, err_vector, err_vector));
        end
        if(params.debug > 1)
            timeinfo = [timeinfo sprintf('check stop:%.2e,', toc - checkstoppartlysmooth)];
        end
    end
    x1 = x2;
    fpre = f1;
    f1 = f2;
    gradf1 = gradf2;
    ng = ng + lng;
    nf = nf + lnf;
    nR = nR + lnR;
    nH = nH + 2;
    times = times + 1;
    X{times + 1} = x1;
    if(times + 1 - params.keepnumX > 0)
        X{times + 1 - params.keepnumX} = [];
    end
    F(times + 1) = f1;
    G(times + 1) = ngf;
    T(times + 1) = toc;
    err_g = G(times + 1);
    if(params.StopCriterion == 0)
        err = 1 - fns.IsStopped(X, F, G, T, params, fns);
    elseif(params.StopCriterion == 1)
        err = err_r;
    elseif(params.StopCriterion == 2)
        err = err_g;
    elseif(params.StopCriterion == 3)
        err = err_g / G(1);
    elseif(params.StopCriterion == 4)
        err = err_v;
    end
    
    if(params.debug > 1)
        cprintf('k', ['\t' timeinfo '\n']);
    end
    if(params.debug > 0)
        if(params.StopCriterion ~= 4)
            fprintf('i:%d,f:%.3e,df/f:%.3e,|gf|:%.3e,initstepsize:%.3e,stepsize:%.3e,time:%.2e,stopls:%d,beta:%.2e,g(s,y):%.2e,g(s,s):%.2e,g(y,y):%.2e\n', ...
                times, f2, err_r, err_g, initstepsize, t, T(end), ls, fns.beta(x1, eta, x2), inp_s_y, inp_s_s, fns.inpro(x2, y, y));
        else
            fprintf('i:%d,f:%.3e,df/f:%.3e,|gf|:%.3e,dx:%.3e,dg:%.3e,initstepsize:%.3e,stepsize:%.3e,time:%.2e,stopls:%d,beta:%.2e,g(s,y):%.2e,g(s,s):%.2e,g(y,y):%.2e\n', ...
                times, f2, err_r, err_g, err_x, err_v, initstepsize, t, T(end), ls, fns.beta(x1, eta, x2), inp_s_y, inp_s_s, fns.inpro(x2, y, y));
        end
    end
end
timecost = toc;
fprintf('Num. of Iter.: %d, time cost: %f seconds, Num. of Fun, Gra, Hv, VT and R: %d, %d, %d, %d, %d \n', times, timecost, nf, ng, nH, nV, nR)
fprintf('f: %e, |gradf|: %e, |gradf/gradf0|: %e, err_v: %e \n', F(end), G(end), G(end)/G(1), err_v)
if(params.isplot)
    figure(1);clf
    scatter(1 : length(F), F, '.');
    ylabel('f(x_i)');
    xlabel('iteration number');
    figure(2);clf
    semilogy(1 : length(G), G, '.');
    ylabel('|grad(x_i)|');
    xlabel('iteration number');
end
varargout{1} = X;
varargout{2} = F;
varargout{3} = G;
varargout{4} = timecost;
varargout{5} = nf;
varargout{6} = ng;
varargout{7} = nH;
varargout{8} = T;
varargout{9} = nV;
varargout{10} = nR;
varargout{11} = times; %------------
end
