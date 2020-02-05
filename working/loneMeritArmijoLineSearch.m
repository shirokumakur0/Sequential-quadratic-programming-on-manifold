function [stepsize, newx] = loneMeritArmijoLineSearch(problem0, rho, xCur, deltaXast, f0, df0, options)

    stepsize = 1;
    newx = problem0.M.retr(xCur, stepsize * deltaXast, rho);
    newf = loneMeritFunction(problem0, newx, rho);
    inner = df0;
    r = 0;

    while newf > f0 - inner  && r <= options.ls_max_steps
        r = r + 1;
        stepsize = stepsize * options.beta;
        inner = inner * stepsize;
        newx = problem0.M.retr(xCur, stepsize * deltaXast, rho);
        newf = loneMeritFunction(problem0, newx, rho);
    end
end
