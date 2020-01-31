function [stepsize, newx] = loneMeritArmijoLineSearch(problem0, rho, xCur, deltaXast, f0, df0, options);

stepsize = 1;
newx = problem0.M.retr(xCur, stepsize * deltaXast, rho);
newf = loneMeritFunction(problem0, newx, rho);
r = 0;

while newf > f0 + df0 && r <= options.ls_max_steps
    r = r + 1;
    stepsize = stepsize * options.beta;
    newx = problem0.M.retr(xCur, stepsize * deltaXast, rho);
    newf = loneMeritFunction(problem0, newx, rho);
end
end
