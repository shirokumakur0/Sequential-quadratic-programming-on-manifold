function [stepsize, newx] = loneMeritArmijoLineSearch(meritproblem, xCur, deltaXast, f0, df0, options)
    stepsize = 1;
    newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
    newf = meritproblem.cost(newx);
    gammadf0 = options.gamma * df0;
    r = 0;
    descriptCost(meritproblem, xCur, deltaXast);
    while newf > ( f0 - gammadf0) && abs(newf - ( f0 - gammadf0)) > 10^(-4)     
        % && r<= options.ls_max_steps
        r = r + 1;
        stepsize = stepsize * options.beta;
        gammadf0 = stepsize * gammadf0;
        newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
        newf = meritproblem.cost(newx);
    end
end
