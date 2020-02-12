function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR, status, timeinfo] = linesearch_exact(fns, eta, x, fv, gf, initslope, initialstepsize, params, timeinfo)
% exactly linesearch algorithm based on fminunc
% locking condition of vector transport is needed.
%
% INPUT:
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
%     fns.gf(x) : return the gradient of objection function at x.
%     fns.Tranv(x1, d, x2, v) : return a tangent vector on T_x2 M, which is given by vector transport that transport v \in T_x M to T_{R_x(d)} M. Here x2 = R_x(d).
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
%     fns.beta(x, eta) : return \|eta\| / \|TR_eta (eta)\|, where TR is the differentiated retraction.
% eta : search direction
% x : current iterate
% fv : function value at the current iterate
% gf : gradient at the current iterate
% alpha : coefficient in the first Wolfe condition
% beta : coefficient in the second Wolfe condition
% minstepsize: allowed minimum step size
% max_step_size: allowed maximum step size
% debug: '0' means almost silence, '1' means some output. '2' means a lot of output, '3' means more than need to know.
%
% OUTPUT:
% new_eta : step size * search direction
% step_size : desired step size
% xnew : next iterate which is R_x(new_eta)
% fnew : function value at the next iterate
% gfnew : gradient at the next iterate
% lnf : number of funtion evaluations in the linesearch algorithm
% lng : number of graduent evaluations in the linesearch algorithm
% lnV : number of vector transport action in the linesearch algorithm
% lnR : number of retraction action in the linesearch algorithm
% TranvParams : Some parameters that may be used in vector transport
% status : 1 means a desired step size is found
%
% By Wen Huang

    alpha = params.alpha;
    beta = params.beta;
    minstepsize = params.minstepsize;
    max_step_size = params.maxstepsize;
    debug = params.debug;
    lnf = 0;
    lng = 0;
    lnV = 0;
    lnR = 0;

    if(debug > 2)
        debugeta1 = fns.STV(x, initialstepsize + 0.00000001, eta);
        debugeta2 = fns.STV(x, initialstepsize, eta);
        debugxnew = fns.R(x, debugeta2);
        [debuggfnew, debugxnew] = fns.gf(debugxnew);
        [debugs, debugeta2, x, debugxnew] = fns.Tranv(x, debugeta2, debugxnew, eta);
        lnV = lnV + 1;
        debugs = fns.STV(debugxnew, 1 / fns.beta(x, debugeta2, debugxnew), debugs);
        newslope = fns.inpro(debugxnew, debuggfnew, debugs);
        fdslope = 100000000*(fns.f(fns.R(x, debugeta1)) - fns.f(debugxnew));
        lnR = lnR + 2;
        lnf = lnf + 1;
        lng = lng + 1;
        cprintf([1,0,1], sprintf('\t\tCheck Locking Condition: slope by inner product / slope by finite different: %e/%e \n', newslope, fdslope));
    end
    
    scaler = fns.inpro(x, gf, gf);
% scaler = 1;
%     options = optimset('TolX', eps, 'TolFun', eps, 'Display', 'off', 'GradObj', 'on');
    options = optimset('Display', 'off', 'GradObj', 'on');
    fhandle = @(t)f_t(t, x, eta, fns, scaler);
    [step_size,fval,exitflag,output] = fminunc(fhandle, initialstepsize, options);
    num = output.funcCount;
    new_eta = fns.STV(x, step_size, eta);
    xnew = fns.R(x, new_eta);
    lnR = lnR + num;
    [fnew, xnew] = fns.f(xnew);
    [gfnew, xnew] = fns.gf(xnew);
    lnf = lnf + num + 1;
    lng = lng + num + 1;
    lnV = lnV + num;
    status = 1;
end

function [f, g] = f_t(t, x, eta, fns, scaler)
    teta = fns.STV(x, t, eta);
    xnew = fns.R(x, teta);
    [f, xnew] = fns.f(xnew);
    f = f / scaler;
    if(nargout > 1)
        [gfnew, xnew] = fns.gf(xnew);
        if(abs(t) < eps)
            g = fns.inpro(x, gfnew, eta) / scaler;
            return;
        end
        
        [s, teta, x, xnew] = fns.Tranv(x, teta, xnew, eta);
        s = fns.STV(x, 1 / fns.beta(x, teta, xnew), s);
        g = fns.inpro(xnew, gfnew, s) / scaler;
    end
end
