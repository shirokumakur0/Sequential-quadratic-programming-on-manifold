function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR, status, timeinfo, newslope, lnftime, lngtime, lnRtime, lnVtime] = linesearch_Strong_Wolfe(fns, eta, x, fv, gf, initslope, initialstepsize, params, timeinfo)

    alpha = params.alpha;
    beta = params.beta;
    minstepsize = params.minstepsize;
    max_step_size = params.maxstepsize;
    debug = params.debug;

    prestepsize = 0;
    fpre = fv;
    newslopepre = initslope;
    step_size = initialstepsize;
    lnf = 0;
    lnftime = 0;
    lng = 0;
    lngtime = 0;
    lnV = 0;
    lnVtime = 0;
    lnR = 0;
    lnRtime = 0;
    status = 1;
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
        lnf = lnf + 1;
        if(fnew > fv + alpha * step_size * initslope)
            [new_eta, step_size, xnew, fnew, gfnew, znf, zng, znV, znR, status, timeinfo, newslope] = Zoom(fns, params, initslope, x, eta, fv, prestepsize, fpre, newslopepre, step_size, fnew);
            lnf = lnf + znf;
            lng = lng + zng;
            lnV = lnV + znV;
            lnR = lnR + znR;
            return;
        end
        
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
        s = fns.STV(xnew, 1 / fns.beta(x, new_eta, xnew), s);
        newslope = fns.inpro(xnew, gfnew, s);
        if(abs(newslope) <= - beta * initslope)
            return;
        end
        if(newslope >= 0)
            [new_eta, step_size, xnew, fnew, gfnew, znf, zng, znV, znR, status, timeinfo, newslope] = Zoom(fns, params, initslope, x, eta, fv, prestepsize, fpre, newslopepre, step_size, fnew);
            lnf = lnf + znf;
            lng = lng + zng;
            lnV = lnV + znV;
            lnR = lnR + znR;
            return;
        end
        prestepsize = step_size;
        fpre = fnew;
        newslopepre = newslope;
        if(step_size >= params.maxstepsize)
            status = -3;
            return;
        end
        step_size = 2 * step_size;
        if(step_size >= params.maxstepsize)
            step_size = params.maxstepsize;
        end
    end
end

function [new_eta, step_size, xnew, fnew, gfnew, znf, zng, znV, znR, status, timeinfo, newslope] = Zoom(fns, params, initslope, x, eta, fv, x1, fx1, slopex1, x2, fx2)
    xlo = x1;
    xhi = x2;
    fxlo = fx1;
    fxhi = fx2;
    znf = 0;
    zng = 0;
    znV = 0;
    znR = 0;
    lnftime = 0;
    lngtime = 0;
    lnVtime = 0;
    lnRtime = 0;
    xloslope = slopex1;
    timeinfo = [];
    status = 1;
    times = 0;
    while (1)
        xdiff = xhi - xlo;
        if(abs(fxhi - fxlo) <= eps || times > 99)
            if(fnew <= fv + params.alpha * step_size * initslope && abs(newslope) <= - 1e-1 * initslope)
                return;
            else
                fnew <= fv + params.alpha * step_size * initslope
                abs(newslope) <= - 1e-1 * initslope
                'h1'
                abs(fxhi - fxlo)
                times
                status = -4;
                return;
            end
        end
        if(abs((fxhi - (fxlo + xloslope * xdiff))) <= eps)
            if(fnew <= fv + params.alpha * step_size * initslope && abs(newslope) <= - 1e-1 * initslope)
                return;
            else
                'h2'
                status = -4;
                return;
            end
        end
        xincr = - xloslope * xdiff * xdiff / 2 / (fxhi - (fxlo + xloslope * xdiff));
        step_size = xlo + xincr;
        new_eta = fns.STV(x, step_size, eta);
        if(params.debug > 1)
            lnRtimestart = toc;
        end
        xnew = fns.R(x, new_eta);
        if(params.debug > 1)
            lnRtime = lnRtime + toc - lnRtimestart;
        end
        znR = znR + 1;
        if(params.debug > 1)
            lnftimestart = toc;
        end
        [fnew, xnew] = fns.f(xnew);
        if(params.debug > 1)
            lnftime = lnftime + toc - lnftimestart;
        end
        znf = znf + 1;
        if(fnew > fv + params.alpha * step_size * initslope)
            xhi = step_size;
            fxhi = fnew;
        else
            if(params.debug > 1)
                lngtimestart = toc;
            end
            [gfnew, xnew] = fns.gf(xnew);
            if(params.debug > 1)
                lngtime = lngtime + toc - lngtimestart;
            end
            zng = zng + 1;
            if(params.debug > 1)
                lnVtimestart = toc;
            end
            [s, new_eta, x, xnew] = fns.Tranv(x, new_eta, xnew, eta);
            if(params.debug > 1)
                lnVtime = lnVtime + toc - lnVtimestart;
            end
            znV = znV + 1;
            s = fns.STV(xnew, 1 / fns.beta(x, new_eta, xnew), s);
            newslope = fns.inpro(xnew, gfnew, s);
            if(abs(newslope) <= - params.beta * initslope)
                status = 1;
                return;
            else
                if(newslope * (xhi - xlo) >= 0)
                    xhi = xlo;
                    fxhi = fxlo;
                end
                xlo = step_size;
                fxlo = fnew;
                xloslope = newslope;
            end
            if(step_size <= params.minstepsize)
                status = -2;
                return;
            end
        end
        times = times + 1;
    end
end
