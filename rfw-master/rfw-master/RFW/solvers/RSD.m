function varargout = RSD(fns,params)
% Solver of Riemannian steepest descent method
%
% INPUT:
% fns : a struct that contains required function handles
% required:
%     fns.f(x) : return objective function value at x.
%     fns.gf(x) : return the gradient of objection function at x.
% not required if params.manifold, params.retraction and params.vector_transport is provided.
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
%     fns.STV(a, v) : return a tangent vector which is a times v.
%     fns.TVaddTV(v1, v2) : return a tangent vector which is equal to v1 adds v2.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
%     fns.beta(x, eta) : return \|eta\| / \|TR_eta (eta)\|, where TR is the differentiated retraction.
%     fns.Tranv(x1, d, x2, v) : return a tangent vector on T_x2 M, which is given by vector transport that transport v \in T_x M to T_{R_x(d)} M. Here x2 = R_x(d).
%
% params : a struct that contains parameters that are used in the solver.
% required:
%     params.x0 : initial approximation of minimizer.
% not required:
%     params.error [1e-5]       : tolerance of stopping criterion.
%     params.StopCriterion [3]  : stopping criterion, 0 means call fns.IsStopped to check
%                                 1 means stop when relative error of objective function is less than tolerance,
%                                 2 means stop when norm of gradient is less than tolerance.
%                                 3 means stop when norm of gradient over initial norm of gradient is less than tolerance,
%     params.max_t [300]        : the maximum number of iterations.
%     params.debug [1]          : '0' means almost silence, '1' means some output. '2' means a lot of output, '3' means more than need to know.
%     params.alpha [1e-4]       : the coefficient of Wolfe first condition or the Armijo condition. It is between (0, 0.5].
%     params.beta [0.999]       : the coefficient of Wolfe second condition. It is between [0.5, 1).
%     params.ratio [0.25]       : the coefficient of the Armijo condition. It is between (0, 1).
%     params.minstepsize [1e-10] : the min step size of line search algorithm. It is between [0, 0.1].
%     params.maxstepsize [200]  : the max step size of line search algorithm. It is between [2, Inf].
%     params.linesearch[1]      : 1 means using inexact line search algorithm with the Armijo conditions and
%                                 2 means using inexact line search algorithm with the Wolfe conditions
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
%     varargout{7} : accumulated computational time for iterations
%     varargout{8} : number of vector transport actions
%     varargout{9} : number of retraction actions
%
% By Wen Huang

fprintf('RSD \n');
if nargin < 2,
    error('Invalid arguments: the number of arguments is invalid.');
end
% set default parameters and check the legal range of parameters
params = set_default_para(params, 'StopCriterion', 3, 'int', 0, 3);
params = set_default_para(params, 'error', 1e-6, 'float', 0, 1);
params = set_default_para(params, 'max_t', 300, 'int', 1, inf);
params = set_default_para(params, 'debug', 1, 'int', 0, 3);
params = set_default_para(params, 'alpha', 1e-4, 'float', 0, 0.5);
params = set_default_para(params, 'beta', 0.999, 'float', 0.5, 1);
params = set_default_para(params, 'minstepsize', 1e-10, 'float', 0, 0.1);
params = set_default_para(params, 'maxstepsize', 200, 'int', 2, inf);
params = set_default_para(params, 'linesearch', 1, 'int', 1, 3);
params = set_default_para(params, 'ratio', 0.25, 'float', 0, 1);
params = set_default_para(params, 'keepnumX', 2100000000, 'int', 1, inf);
params = set_default_para(params, 'isplot', 0, 'int', 0, 1);
params = set_default_para(params, 'initstepsize', 1, 'float', 0, inf);
if(isfield(params, 'manifold') && isfield(params, 'retraction') && isfield(params, 'vector_transport'))
    fs = choose_manifold(params.manifold, params.retraction, params.vector_transport);
    fns.R = fs.R;
    fns.STV = fs.STV;
    fns.TVaddTV = fs.TVaddTV;
    fns.beta = fs.beta;
    fns.inpro = fs.inpro;
    fns.Tranv = fs.Tranv;
end
% Check whether function handles are legal
check_function_handle(fns, 'f');
check_function_handle(fns, 'gf');
check_function_handle(fns, 'R');
check_function_handle(fns, 'STV');
check_function_handle(fns, 'TVaddTV');
check_function_handle(fns, 'beta');
check_function_handle(fns, 'inpro');
check_function_handle(fns, 'Tranv');
if(params.StopCriterion == 0)
    check_function_handle(fns, 'IsStopped');
end
% Check whether required parameters are legal
if(~isfield(params, 'x0'))
    error('Invalid arguments: missing initial x0');
end
% Initialization


err = inf;
times = 0;
x1 = params.x0;
[f1, x1] = fns.f(x1);
[gradf1, x1] = fns.gf(x1);

X{1} = x1;
F(1) = f1;
G(1) = sqrt(fns.inpro(x1, gradf1, gradf1));
err_g = inf;

nf = 1;
ng = 1;
nV = 0;
nR = 0;
fprintf('times : %d, f(x) : %e \n', times, f1);
timeinfo = '';
tic
T(1) = 0;
Tf(1) = 0;
Tg(1) = 0;
Tr(1) = 0;
Tv(1) = 0;
while(err > params.error && times < params.max_t)
    eta = fns.STV(x1, -1, gradf1);
    if(params.debug > 1)
        timeinfo = '';
        linesearchtimestart = toc;
    end
    initslope = fns.inpro(x1, eta, gradf1);
    if(times == 0)
        initstepsize = params.initstepsize;
    else
       initstepsize = 2 * (f1 - fpre) / initslope; 
    end
    if(params.linesearch == 1)
        if(err_g < params.accuracy)
            [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, lnftime, lngtime, lnRtime, lnVtime] = choose_step_size(fns, eta, x1, f1, gradf1, params);
        else
            [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, timeinfo, lnftime, lngtime, lnRtime, lnVtime] = linesearch_Armijo(fns, eta, x1, f1, gradf1, initslope, initstepsize, params, timeinfo);
        end
    elseif(params.linesearch == 2)
        if(err_g < params.accuracy)
            [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR,ls, lnftime, lngtime, lnRtime, lnVtime] = choose_step_size(fns, eta, x1, f1, gradf1, params);
        else
            [eta, t, x2, f2, gradf2, lnf, lng, lnV, lnR, ls, timeinfo, lnftime, lngtime, lnRtime, lnVtime] = linesearch_Wolfe(fns, eta, x1, f1, gradf1, initslope, initstepsize, params, timeinfo);
        end
    end
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('linesearch total:%.2e,', toc - linesearchtimestart)];
    end
    if(ls < 0)
        fprintf('warning: line search fails!\n');
        break; 
    end
    % check stopping cretirion
    err_r = abs((f2 - f1) / f2);
    x1 = x2;
    fpre = f1;
    f1 = f2;
    gradf1 = gradf2;
    ng = ng + lng;
    nf = nf + lnf;
    nV = nV + lnV;
    nR = nR + lnR;
    times = times + 1;
    X{times + 1} = x1;
    if(times + 1 - params.keepnumX > 0)
        X{times + 1 - params.keepnumX} = [];
    end
    F(times + 1) = f1;
    G(times + 1) = sqrt(fns.inpro(x1, gradf1, gradf1));
    T(times + 1) = toc;
    err_g = G(times + 1);
    Tf(times + 1) = Tf(times) + lnftime;
    Tg(times + 1) = Tg(times) + lngtime;
    Tr(times + 1) = Tr(times) + lnRtime;
    Tv(times + 1) = Tv(times) + lnVtime;

    if(params.StopCriterion == 0)
        err = 1 - fns.IsStopped(X, F, G, T, params, fns);
    elseif(params.StopCriterion == 1)
        err = err_r;
    elseif(params.StopCriterion == 2)
        err = err_g;
    elseif(params.StopCriterion == 3)
        err = err_g / G(1);
    end
    if(params.debug > 1)
        cprintf('k', ['\t' timeinfo '\n']);
    end
    if(params.debug > 0)
        fprintf('i:%d,f:%.3e,df/f:%.3e,|gf|:%.3e,initstepsize:%.3e,stepsize:%.3e,stopls:%d,time:%.2e\n', ...
            times, f2, err_r, err_g, initstepsize, t, ls, T(end));
    end
end
timecost = toc;
fprintf('Num. of Iter.: %d, time cost: %f seconds, Num. of Fun, Gra, VT and R: %d, %d, %d, %d \n', times, timecost, nf, ng, nV, nR)
fprintf('f: %e, |gradf|: %e, |gradf/gradf0|: %e \n', F(end), G(end), G(end)/G(1))
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
varargout{7} = T;
varargout{8} = nV;
varargout{9} = nR;
varargout{10} = times;
varargout{11} = Tf;
varargout{12} = Tg;
varargout{13} = Tr;
varargout{14} = Tv;
end
