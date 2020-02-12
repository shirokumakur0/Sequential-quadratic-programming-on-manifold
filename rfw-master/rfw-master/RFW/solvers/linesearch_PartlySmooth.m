function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR, status, timeinfo] = linesearch_PartlySmooth(fns, eta, x, fv, gf, initslope, initialstepsize, params, timeinfo)
% This function apply linesearch algorithm to find a point for partly smooth functions
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
% status : -4 means the maximum number of iteration is reached, 
%          -3 means step size excesses the maximum step size
%          -2 means step size reaches the minimum step size
%          -1 means unable to find a point that satisfies the curvature condition
%           1 means a desired step size is found
%
% By Wen Huang

    alpha = params.alpha;
    beta = params.beta;
    minstepsize = params.minstepsize;
    max_step_size = params.maxstepsize;
    debug = params.debug;
    lnf = 0;
    lnftime = 0;
    lng = 0;
    lngtime = 0;
    lnV = 0;
    lnVtime = 0;
    lnR = 0;
    lnRtime = 0;
%     initslope = fns.inpro(x, eta, gf);
    a = 0;
    b = inf;
    step_size = initialstepsize;
    times = 0;
    while(1)
        new_eta = fns.STV(x, step_size, eta);
        if(params.debug > 1)
            lnRtimestart = toc;
        end
        xnew = fns.R(x, new_eta);
        if(params.debug > 1)
            lnRtime = lnRtime + toc - lnRtimestart;
        end
        lnR = lnR + 1;
        if(params.debug > 1)
            lnftimestart = toc;
        end
        [fnew, xnew] = fns.f(xnew);
        if(params.debug > 1)
            lnftime = lnftime + toc - lnftimestart;
        end
        lnf = lnf + 1;
        if(params.debug > 1)
            lngtimestart = toc;
        end
        [gfnew, xnew] = fns.gf(xnew);
        if(params.debug > 1)
            lngtime = lngtime + toc - lngtimestart;
        end
        lng = lng + 1;
        if(params.debug > 1)
            lnVtimestart = toc;
        end
        [s, new_eta, x, xnew] = fns.Tranv(x, new_eta, xnew, eta);
        if(params.debug > 1)
            lnVtime = lnVtime + toc - lnVtimestart;
        end
        lnV = lnV + 1;
        s = fns.STV(x, 1 / fns.beta(x, new_eta, xnew), s);
        newslope = fns.inpro(xnew, gfnew, s);
        if(debug > 2)
            debugeta1 = fns.STV(x, step_size + 0.00000001, eta);
            debugeta2 = fns.STV(x, step_size, eta);
            if(params.debug > 1)
                lnRtimestart = toc;
            end
            fdslope = 100000000*(fns.f(fns.R(x, debugeta1)) - fns.f(fns.R(x, debugeta2)));
            if(params.debug > 1)
                lnRtime = lnRtime + toc - lnRtimestart;
            end
            lnR = lnR + 2;
            cprintf([1,0,1],sprintf('\t\tAnalytical gradient: %e, finite different gradient: %e, they should be similar \n', newslope, fdslope));
        end
        if(fnew - fv >= alpha * initslope * step_size)
            b = step_size;
        elseif(newslope <= beta * initslope || ~fns.gf_exists(xnew))
            a = step_size;
        else
            new_eta = fns.STV(x, step_size, eta);
            status = 1;
            break;
        end
        if(b < inf)
            step_size = (a + b) / 2;
        else
            step_size = 2 * a;
        end
        if(step_size < minstepsize)
            fprintf('warning: linesearch reaches the min step size, it may need to be diminished \n')
            new_eta = fns.STV(x, step_size, eta);
            status = -2;
            if(params.debug > 1)
                timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
            end
            return;
        elseif(step_size > max_step_size)
            fprintf('warning: linesearch reaches the max step size, it may need to be enlarged \n')
            new_eta = fns.STV(x, step_size, eta);
            status = -3;
            if(params.debug > 1)
                timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
            end
            return;
        end
        times = times + 1;
        if(times > 200)
            fprintf('warning: linesearch reaches the max num of iters \n')
            new_eta = fns.STV(x, step_size, eta);
            status = -4;
            break;
        end
    end
    if(params.debug > 1)
        timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
    end
end
