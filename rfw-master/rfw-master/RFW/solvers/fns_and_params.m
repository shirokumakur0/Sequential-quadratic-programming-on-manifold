
fns : a struct that contains required function handles
    fns.f(x) : return objective function value at x.
    fns.gf(x) : return the gradient of objection function at x.
    fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
    fns.Hv(x, H, v) : return tangent vector which is given by that operator H affect on tangent vector v on T_x M.
    fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
    fns.beta(x, eta) : return \|eta\| / \|TR_eta (eta)\|, where TR is the differentiated retraction.
    fns.phi(x, y, s, B, inp_s_y, sBs, fns, v, times) : return the coefficient of RBroyden Family update formula.
    fns.Tranv(x1, d, x2, v) : return a tangent vector on T_x2 M, which is given by vector transport that transport v \in T_x M to T_{R_x(d)} M. Here x2 = R_x(d).
    fns.TranH(x1, d, x2, H1) : return a operator, H2, that affects on T_x2 M. H2 satisfies that for any v \in T_x2(M), H2(v) = Tran (H1 (Tran^{-1}(v)))
                            where Tran is the vector transport fns.Tranv(x1, d, x2, v).
    fns.rank1operator(x, v1, v2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = g_x(v2, v) * v1.
    fns.operadd(x, H1, H2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = H1(v) + H2(v).
    fns.Tranv_Params(x1, d, x2, v, TranvParams) : return the tangent vector which is the same as fns.Tranv(x1, d, x2, v), 
                                                  TranvParams is some parameters such that the cost becomes smaller and 
                                                  are given by fns.Tranv or fns.invTranv.
    fns.proj(x, eta) : return P_x(eta) by projecting v to tangent space of x.
    fns.coTVtimesTV(x, coeta, eta) : return a number which is the result of the cotangent vector, coeta, acts on the tangent vector, eta.
    fns.coTV(x1, d, x2, v) : return the desired cotangent vector required by Ring and Wirth's RBFGS
    fns.SRTV(x, m, seed) : return m random tangent vectors in T_x M, seed is the random seed.

params : a struct that contains parameters that are used in the solver.
    params.x0 : initial approximation of minimizer.
    params.error [1e-5]       : tolerance of stopping criterion.
    params.StopCriterion [2]  : stopping criterion, 1 means stop when relative error of objective function is less than tolerance,
                                2 means stop when norm of gradient is less than tolerance,
                                3 means stop when norm of gradient over initial norm of gradient is less than tolerance,
                                4 means stopping criterion for partly smooth function.
    params.m [4]              : number of s and y for storage
    params.restart [1]        : 1 means restart if curvature condition is not satisfied, 2 means return directly.
                                (This should not happen theoretically. In practice, only numerical error make it occur occasionally.)
    params.max_t [300]        : the maximum number of iterations.
    params.debug [1]          : '0' means almost silence, '1' means some output. '2' means a lot of output, '3' means more than need to know.
    params.alpha [1e-4]       : the coefficient of Wolfe first condition. It is between (0, 0.5].
    params.beta [0.999]       : the coefficient of Wolfe second condition. It is between [0.5, 1).
    params.minstepsize [1e-4] : the min step size of line search algorithm. It is between [0, 0.1].
    params.maxstepsize [200]  : the max step size of line search algorithm. It is between [2, Inf].
    params.err_x [1e-5]       : tolerance for checking whether iterates are close to convergence. It is between [0, 1].
    params.num_Grads [10]     : number of gradients for checking whether iterates are close to Clark stationary point. It is between [0, inf].
    params.theta[0.1]         : A stopping criterion parameter for solving local model
    params.kappa[0.9]         : A stopping criterion parameter for solving local model
    params.min_innit[0]       : the minimum number of iterations for solving local model
    params.max_innit[200]     : the maxmum number of iterations for solving local model
    params.linesearch[1]      : 1 means using exact line search algorithm and 2 means using inexact line search algorithm with the Wolfe conditions
    params.linesearch[1]      : 1 means using inexact line search algorithm with the Wolfe conditions and 
                                2 means using inexact line search algorithm with the Armijo conditions 
FOR GS
    params.nu_err[1e-10]      : tolerance of optimality tolerance
    params.nd_err[1e-10]      : tolerance for checking whether it is close to a Clark stationary point
    params.ns [20]             : the number of samples
    params.epsilon [1e-10]    : sampling radius
    params.nu [0]             : Initial optimality tolerance
    params.theta [1]          : Optimality tolerance reduction factor
    params.mu [1]             : Sampling radius reduction factor
    params.gsbeta [0]         : Armijo parameter
    params.gamma [0.5]        : Backtracking reduction factor
    params.max_ls [40]        : maximum number of iteration for backtracking in linesearch


The following give some examples of manifold. For each manifold, some retractions and vector transports are provided.
    params.manifold : the manifold on which objective function is defined.
        '1', sphere : S^{n - 1}
        '2', Stiefel manifold : St(p, n)
        '3', orthogonal group : O(n)
        '4', grassmann manifold : Gr(p, n)

varargout : contains output parameters
    varargout{1} : all iterates
    varargout{2} : function values of all iterates
    varargout{3} : norm of gradient of all iterates
    varargout{4} : total computational time
    varargout{5} : number of function evaluations
    varargout{6} : number of gradient evaluations
    varargout{7} : accumulated computational time for iterations
    varargout{8} : number of vector transport actions without parameters
    varargout{9} : number of vector transport action with parameters
    varargout{10} : number of retraction actions

varargout : contains output parameters
    varargout{1} : all iterates
    varargout{2} : function values of all iterates
    varargout{3} : norm of gradient of all iterates
    varargout{4} : total computational time
    varargout{5} : number of function evaluations
    varargout{6} : number of gradient evaluations
    varargout{7} : number of Hessian actions
    varargout{8} : accumulated computational time for iterations
    varargout{9} : number of vector transport actions
    varargout{10} : number of retraction actions

the solver of local model is a modified version of that in http://www.math.fsu.edu/~cbaker/GenRTR/




        '1', noniso: vector transport by projection
        '2', Emb Iso: parallel transport
        '3', Emb Iso: direct-rotation with basis on tangent space
        '4', Emb Iso: direct-rotation with basis on normal space
        '5', Emb Iso: smooth basis on tangent space
        '6', Emb Iso: smooth basis on normal space
        '7', Emb Iso: intrinsic method
        '8', Quo Iso: parallel transport
        '9', Quo Iso: direct-rotation with basis on horizontal space
        '10', Quo Iso: direct-rotation with basis on vertical space
        '11',Quo Iso: smooth basis on horizontal space
        '12',Quo Iso: smooth basis on vertical space
        '13',Quo Iso: intrinsic method
        '14',Emb Locking-Iso(qf): direct-rotation with basis on tangent space
        '15',Emb Locking-Iso(qf): direct-rotation with basis on normal space
        '16',Emb Locking-Iso(qf): smooth basis on tangent space
        '17',Emb Locking-Iso(qf): smooth basis on normal space
        '18',Emb Locking-Iso(qf): construct vector transport
        '19',Emb Locking-Iso(qf): intrinsic method
        '20',Quo Locking-Iso(qf): direct-rotation with basis on horizontal space
        '21',Quo Locking-Iso(qf): direct-rotation with basis on vertical space
        '22',Quo Locking-Iso(qf): smooth basis on horizontal space
        '23',Quo Locking-Iso(qf): smooth basis on vertical space
        '24',Quo Locking-Iso(qf): construct vector transport
        '25',Quo Locking-Iso(qf): intrinsic method
        '26',Emb Locking-nonIso(qf): vector transport by rojection
