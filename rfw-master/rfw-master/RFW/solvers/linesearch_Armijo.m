function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR, status, timeinfo, lnftime, lngtime, lnRtime, lnVtime] = linesearch_Armijo(fns, eta, x, fv, gf, initslope, initialstepsize, params, timeinfo)
% This function apply linesearch algorithm to find a point which satisfies the Amijo condition.
%
% INPUT:
% fns : a struct that contains required function handles
%     fns.f(x) : return objective function value at x.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
%     fns.gf(x) : return the gradient of objection function at x.
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
% eta : search direction
% x : current iterate
% fv : function value at the current iterate
% gf : gradient at the current iterate
% alpha : coefficient in the Amijo condition
% ratio : coefficient in the Amijo condition
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
% status : -4 means some serious errors happend,
%          -3 means step size excesses the maximum step size
%          -2 means step size reaches the minimum step size
%          -1 means unable to find a point that satisfies the curvature condition
%           1 means a desired step size is found
%
% By Wen Huang
global issimple;
step_size = initialstepsize;
alpha = params.alpha;
ratio = params.ratio;
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
m = 0;
new_eta = fns.STV(x, step_size, eta);

if(issimple == 1)
    
    lnRtimestart = toc;
    xnew = fns.R(x, new_eta);
        lnRtime = lnRtime + toc - lnRtimestart;
    lnR = lnR + 1;
    while(min(eig(xnew.U)) <=  0)
        step_size = step_size * ratio;
        new_eta = fns.STV(x, step_size, eta);
        
        lnRtimestart = toc;
        xnew = fns.R(x, new_eta);
        lnRtime = lnRtime + toc - lnRtimestart;
        lnR = lnR + 1;
        fprintf('warning: negative eig, pds:%e \n', step_size); %------
    end
end



lnRtimestart = toc;
xnew = fns.R(x, new_eta);
lnRtime = lnRtime + toc - lnRtimestart;
lnR = lnR + 1;

lnftimestart = toc;
fnew = fns.f(xnew);
lnftime = lnftime + toc - lnftimestart;
lnf = lnf + 1;
times = 0;
while(fv - fnew < - alpha * initslope * step_size)
    step_size = ratio * step_size;
    new_eta = fns.STV(x, step_size, eta);
    if(params.debug > 1)
        lnRtimestart = toc;
    end
    
    lnRtimestart = toc;
    xnew = fns.R(x, new_eta);
    lnRtime = lnRtime + toc - lnRtimestart;
    
    if(params.debug > 1)
        lnRtime = lnRtime + toc - lnRtimestart;
    end
    lnR = lnR + 1;
    
    lnftimestart = toc;
    [fnew, xnew] = fns.f(xnew);
    lnftime = lnftime + toc - lnftimestart;
    lnf = lnf + 1;
    
    if(step_size < minstepsize)
        step_size = minstepsize;
        
        lnRtimestart = toc;
        xnew = fns.R(x, fns.STV(x, step_size, eta));
        lnRtime = lnRtime + toc - lnRtimestart;
        lnR = lnR + 1;
        
        lnftimestart = toc;
        [fnew, xnew] = fns.f(xnew);
        lnftime = lnftime + toc - lnftimestart;
        lnf = lnf + 1;
        
        lngtimestart = toc;
        [gfnew, xnew] = fns.gf(xnew);
        lngtime = lngtime + toc - lngtimestart;
        lng = lng + 1;
        
        new_eta = fns.STV(x, minstepsize, eta);
        status = -2;
        if(debug > 0)
            fprintf('warning: linesearch reaches the min step size, it may need to be diminished \n')
        end
        if(params.debug > 1)
            timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
        end
        return;
    end
    times = times + 1;
end
status = 1;

lngtimestart = toc;
[gfnew, xnew] = fns.gf(xnew);
lngtime = lngtime + toc - lngtimestart;
lng = lng + 1;

if(params.debug > 1)
    timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
end
end
