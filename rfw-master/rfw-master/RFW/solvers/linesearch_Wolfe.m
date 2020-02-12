function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR, status, timeinfo, lnftime, lngtime, lnRtime, lnVtime] = linesearch_Wolfe(fns, eta, x, fv, gf, initslope, initialstepsize, params, timeinfo)
% This function apply linesearch algorithm to find a point which satisfies Wolfe conditions.
% See section 6.3.2 & Algorithm 6.3.1mod of J.E. Dennis, Jr. Robert B. Schnabel's book :
% Numerical Methods for Unconstrained Optimization and Nonlinear Equation.
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
% status : -4 means some serious errors happend,
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

step_size = initialstepsize;
lnf = 0;
lnftime = 0;
lng = 0;
lngtime = 0;
lnV = 0;
lnVtime = 0;
lnR = 0;
lnRtime = 0;
while(1)
    if(isnan(step_size) || isinf(step_size))
        status = -4;
        new_eta = 0;
        gfnew = gf;
        if(params.debug > 1)
            timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
        end
        return;
    end
    new_eta = fns.STV(x, step_size, eta);
    
    lnRtimestart = toc;
    xnew = fns.R(x, new_eta);
    lnRtime = lnRtime + toc - lnRtimestart;
    lnR = lnR + 1;
    
    lnftimestart = toc;
    [fnew, xnew] = fns.f(xnew);
    lnftime = lnftime + toc - lnftimestart;
    lnf = lnf + 1;
    
    if(fnew <= fv + alpha * step_size * initslope)
        
        lngtimestart = toc;
        [gfnew, xnew] = fns.gf(xnew);
        lngtime = lngtime + toc - lngtimestart;
        new_eta = fns.STV(x, step_size, eta);
        lng = lng + 1;
        
        lnVtimestart = toc;
        [s, new_eta, x, xnew] = fns.Tranv(x, new_eta, xnew, eta);
        lnVtime = lnVtime + toc - lnVtimestart;
        lnV = lnV + 1;
        s = fns.STV(xnew, 1 / fns.beta(x, new_eta, xnew), s);
        newslope = fns.inpro(xnew, gfnew, s);
        if(debug > 2)
            debugeta1 = fns.STV(x, step_size + 0.00000001, eta);
            debugeta2 = fns.STV(x, step_size, eta);
            
            lnRtimestart = toc;
            debugx1 = fns.R(x, debugeta1);
            debugx2 = fns.R(x, debugeta2);
            lnRtime = lnRtime + toc - lnRtimestart;
            
            lnftimestart = toc;
            fdslope = 100000000*(fns.f(debugx1) - fns.f(debugx2));
            lnftime = lnftime + toc - lnftimestart;
            
            lnR = lnR + 2;
            lnf = lnf + 2;
            cprintf([1,0,1], sprintf('\t\tCheck Locking Condition: slope by inner product / slope by finite different: %e/%e \n', newslope, fdslope));
        end
        if(newslope < beta * initslope)
            times = 0;
            while(fnew <= fv + alpha * step_size * initslope && newslope < beta * initslope && step_size < max_step_size)
                pre_step_size = step_size;
                fnew_pre = fnew;
                step_size = min(2 * step_size, max_step_size);
                new_eta = fns.STV(x, step_size, eta);
                
                lnRtimestart = toc;
                xnew = fns.R(x, new_eta);
                lnRtime = lnRtime + toc - lnRtimestart;
                lnR = lnR + 1;
                
                lnftimestart = toc;
                [fnew, xnew] = fns.f(xnew);
                lnftime = lnftime + toc - lnftimestart;
                lnf = lnf + 1;
                
                if(fnew <= fv + alpha * step_size * initslope)
                    
                    lngtimestart = toc;
                    [gfnew, xnew] = fns.gf(xnew);
                    lngtime = lngtime + toc - lngtimestart;
                    lng = lng + 1;
                    
                    lnVtimestart = toc;
                    [s, new_eta, x, xnew] = fns.Tranv(x, new_eta, xnew, eta);
                    lnVtime = lnVtime + toc - lnVtimestart;
                    lnV = lnV + 1;
                    
                    s = fns.STV(xnew, 1 / fns.beta(x, new_eta, xnew), s);
                    newslope = fns.inpro(xnew, gfnew, s);
                end
                times = times + 1;
                if(times > 10)
                    break;
                end
            end
            if(step_size == max_step_size)
                
                lngtimestart = toc;
                [gfnew, xnew] = fns.gf(xnew);
                lngtime = lngtime + toc - lngtimestart;
                lng = lng + 1;
                new_eta = fns.STV(x, step_size, eta);
                status = -3;
                if(debug > 0)
                    fprintf('warning: linesearch reaches the max step size, it may need to be enlarged \n')
                end
                if(params.debug > 1)
                    timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
                end
                return;
            end
            if(step_size ~= initialstepsize && fnew > fv + alpha * step_size * initslope)
                step_sizelo = min(step_size, pre_step_size);
                step_sizediff = abs(pre_step_size - step_size);
                if(step_size < pre_step_size)
                    flo = fnew;
                    fhi = fnew_pre;
                else
                    flo = fnew_pre;
                    fhi = fnew;
                end
                times = 0;
                while(newslope < beta * initslope && step_sizediff >= minstepsize)
                    step_sizeincr = - newslope * step_sizediff^2 / 2 / (fhi - (flo + newslope * step_sizediff));
                    if(step_sizeincr < 0.2 * step_sizediff)
                        step_sizeincr = 0.2 * step_sizediff;
                    end
                    step_size = step_sizelo + step_sizeincr;
                    new_eta = fns.STV(x, step_size, eta);
                    
                    lnRtimestart = toc;
                    xnew = fns.R(x, new_eta);
                    lnRtime = lnRtime + toc - lnRtimestart;
                    lnR = lnR + 1;
                    
                    lnftimestart = toc;
                    [fnew, xnew] = fns.f(xnew);
                    lnftime = lnftime + toc - lnftimestart;
                    lnf = lnf + 1;
                    
                    if(fnew > fv + alpha * step_size * initslope)
                        step_sizediff = step_sizeincr;
                        fhi = fnew;
                    else
                        
                        lngtimestart = toc;
                        [gfnew, xnew] = fns.gf(xnew);
                        lngtime = lngtime + toc - lngtimestart;
                        lng = lng + 1;
                        
                        lnVtimestart = toc;
                        [s, new_eta, x, xnew] = fns.Tranv(x, new_eta, xnew, eta);
                        lnVtime = lnVtime + toc - lnVtimestart;
                        lnV = lnV + 1;
                        
                        s = fns.STV(xnew, 1 / fns.beta(x, new_eta, xnew), s);
                        newslope = fns.inpro(xnew, gfnew, s);
                        if(newslope < beta * initslope)
                            step_sizelo = step_size;
                            step_sizediff = step_sizediff - step_sizeincr;
                            flo = fnew;
                        end
                    end
                    times = times + 1;
                    if(times > 10)
                        break;
                    end
                end
                if(newslope < beta * initslope)
                    fnew = flo;
                    new_eta = fns.STV(x, step_sizelo, eta);
                    
                    lnRtimestart = toc;
                    xnew = fns.R(x, new_eta);
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
                    
                    step_size = step_sizelo;
                    status = -1;
                    if(debug > 0)
                        fprintf('warning: linesearch couldn''t find a point that satisfies curvature condition, \n\t\t newslope: %f, initslope: %f, newslope/initslope: %f\n', newslope, initslope, newslope/initslope)
                    end
                    if(params.debug > 1)
                        timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
                    end
                    return;
                end
            end
        end
        
        status = 1;
        if(params.debug > 1)
            timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
        end
        return;
    end
    if(step_size < minstepsize)
        step_size = minstepsize;
        new_eta = fns.STV(x, step_size, eta);
        
        lnRtimestart = toc;
        xnew = fns.R(x, new_eta);
        lnRtime = lnRtime + toc - lnRtimestart;
        lnR = lnR + 1;
        
        lnftimestart = toc;
        [fnew, xnew] = fns.f(xnew);
        lnftime = lnftime + toc - lnftimestart;
        
        lngtimestart = toc;
        [gfnew, xnew] = fns.gf(xnew);
        lngtime = lngtime + toc - lngtimestart;
        lnf = lnf + 1;
        lng = lng + 1;
        %             new_eta = minstepsize * eta;
        status = -2;
        if(debug > 0)
            fprintf('warning: linesearch reaches the min step size, it may need to be diminished \n')
        end
        if(params.debug > 1)
            timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
        end
        return;
    else
        if(step_size == initialstepsize)
            steptemp = - initslope / 2 / (fnew - fv - initslope);
        else
            n1 = 1 / step_size / step_size;
            n2 = 1 / pre_step_size / pre_step_size;
            M = [n1, -n2; - pre_step_size * n1, step_size * n2] * [fnew - fv - step_size * initslope; fnew_pre - fv - pre_step_size * initslope] / (step_size - pre_step_size);
            a = M(1);
            b = M(2);
            disc = b^2 - 3 * a * initslope;
            if(abs(a) < 1e-10)
                steptemp = -initslope / 2 / b;
            else
                steptemp = (-b + sqrt(disc)) / 3 / a;
            end
            if(steptemp > 0.5 * step_size)
                steptemp = 0.5 * step_size;
            end
        end
        pre_step_size = step_size;
        fnew_pre = fnew;
        if(steptemp <= 0.1 * step_size)
            step_size = 0.1 * step_size;
        else
            step_size = steptemp;
        end
    end
end
if(params.debug > 1)
    timeinfo = [timeinfo sprintf('[linesearch:lnftime(%d:%.2e),lngtime(%d:%.2e),lnVtime(%d:%.2e),lnRtime(%d:%.2e)]', lnf, lnftime, lng, lngtime, lnV, lnVtime, lnR, lnRtime)];
end
end