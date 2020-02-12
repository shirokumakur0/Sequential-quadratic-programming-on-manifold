function varargout = check_grad_Hess(fns,params)
% Solver of Riemannian trust region steepest descent method
%
% INPUT: 
% fns : a struct that contains required function handles
% required:
%     fns.f(x) : return objective function value at x.
%     fns.gf(x) : return the gradient of objection function at x.
% not required if params.manifold, params.retraction and params.vector_transport is provided.
%     fns.STV(a, v) : return a tangent vector which is a times v.
%     fns.TVaddTV(v1, v2) : return a tangent vector which is equal to v1 adds v2.
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
%     fns.proj(x, eta) : return P_x(eta) by projecting v to tangent space of x.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
% 
% params : a struct that contains parameters that are used in the solver.
% required:
%     params.x0 : initial approximation of minimizer.
% not required:
%     params.error [1e-5]       : tolerance of stopping criterion.
%     params.StopCriterion [2]  : stopping criterion, 1 means stop when relative error of objective function is less than tolerance,
%                                 2 means stop when norm of gradient is less than tolerance,
%                                 3 means stop when norm of gradient over initial norm of gradient is less than tolerance,
%                                 4 means stopping criterion for partly smooth function.
%     params.max_t [300]        : the maximum number of iterations.
%     params.debug [1]          : '0' means almost silence, '1' means some output. '2' means a lot of output, '3' means more than need to know.
%     params.c [0.1]            : acceptance or rejection constant
%     params.tau1 [0.25]        : shrink radius constant
%     params.tau2 [2]           : expansion radius constant
%     params.max_Delta [200]    : max of Delta
%     params.min_Delta [1e-10]   : min of Delta
%     params.err_x [1e-5]       : tolerance for checking whether iterates are close to convergence. It is between [0, 1].
%     params.num_Grads [10]     : number of gradients for checking whether iterates are close to Clark stationary point. It is between [0, inf].
%     params.theta[0.1]         : A stopping criterion parameter for solving local model
%     params.kappa[0.9]         : A stopping criterion parameter for solving local model
%     params.min_innit[0]       : the minimum number of iterations for solving local model
%     params.max_innit[200]     : the maxmum number of iterations for solving local model
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
%     varargout{8} : number of retraction actions
% 
% By Wen Huang, the solver of local model is a modified version of that in http://www.math.fsu.edu/~cbaker/GenRTR/

    fprintf('check grad Hess\n')
    if nargin < 2,
       error('Invalid arguments: the number of arguments is invalid.');
    end
    % set default parameters and check the legal range of parameters
    params = set_default_para(params, 'error', 1e-5, 'float', 0, 1);
    params = set_default_para(params, 'err_x', 1e-5, 'float', 0, 1);
    params = set_default_para(params, 'num_Grads', 10, 'int', 0, inf);
    params = set_default_para(params, 'max_t', 300, 'int', 1, inf);
    params = set_default_para(params, 'debug', 1, 'int', 0, 3);
    params = set_default_para(params, 'c', 0.1, 'float', 0, 0.5);
    params = set_default_para(params, 'tau1', 0.25, 'float', 0.1, 0.9);
    params = set_default_para(params, 'tau2', 2, 'float', 1.5, 5);
    params = set_default_para(params, 'max_Delta', 200, 'float', 1, inf);
    params = set_default_para(params, 'min_Delta', 1e-10, 'float', 0, 1);
    params = set_default_para(params, 'theta', 0.1, 'float', 0, 1);
    params = set_default_para(params, 'kappa', 0.9, 'float', 0, 1);
    params = set_default_para(params, 'min_innit', 0, 'int', 0, inf);
    params = set_default_para(params, 'max_innit', 200, 'int', 0, inf);
    if(isfield(params, 'manifold') && isfield(params, 'retraction') && isfield(params, 'vector_transport'))
        fs = choose_manifold(params.manifold, params.retraction, params.vector_transport);
        fns.R = fs.R;
        fns.proj = fs.proj;
        fns.inpro = fs.inpro;
        fns.STV = fs.STV;
        fns.TVaddTV = fs.TVaddTV;
    end
    % Check whether function handles are legal
    check_function_handle(fns, 'f');
    check_function_handle(fns, 'gf');
    check_function_handle(fns, 'hessianv');
    check_function_handle(fns, 'STV');
    check_function_handle(fns, 'TVaddTV');
    check_function_handle(fns, 'R');
    check_function_handle(fns, 'proj');
    check_function_handle(fns, 'inpro');
    % Check whether required parameters are legal
    if(~isfield(params, 'x0'))
        error('Invalid arguments: missing initial x0');
    end

    fns.g = fns.inpro;

    times = 0;
    x = params.x0;
    [fx, x] = fns.f(x);
    [gradf, x] = fns.gf(x);
    RandTV.v = rand(size(gradf.v));
    eta = fns.proj(x, RandTV);    % randomly choose a tangent vector
    neta = fns.g(x, eta, eta);
    eta = fns.STV(x, 100 / sqrt(neta), eta); % initial length of eta is 100.

    % the length of eta variance from 100 to 100*2^(-30) approx 9e-8
    t = 1;
    for i = 1 : 30
        t = t * 0.5;
        teta = fns.STV(x, t, eta);
        xtemp = fns.R(x, teta);
        [fxtemp, xtemp] = fns.f(xtemp);
        YY(i) = log(fxtemp - fx - fns.g(x, gradf, teta) - 0.5 * fns.g(x, teta, fns.hessianv(x, teta)));
        YX(i) = log(sqrt(fns.g(x, teta, teta)));
        fprintf('i:%d, x:%e, %e, %e, %e, ', i, YX(i), fxtemp - fx, fns.g(x, gradf, teta), (fxtemp - fx) / fns.g(x, gradf, teta));
        fprintf('%e, %e, \n', (fxtemp - fx - fns.g(x, gradf, teta)), (fxtemp - fx - fns.g(x, gradf, teta)) / fns.g(x, teta, fns.hessianv(x, teta)));
    end
    plot(real(YX), real(YY));
    
    fprintf('If the retraction is second order or the retraction is not second order but x is at a critical point\n')
    fprintf('then the slope of the curve should be approximately 3.\n');

    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = 0;
    varargout{5} = [];
    varargout{6} = [];
    varargout{7} = [];
    varargout{8} = [];
end
