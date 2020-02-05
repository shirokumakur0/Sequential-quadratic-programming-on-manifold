function [meritproblem, f0, df0] = makeMeritProblem(problem0, xCur, rho, deltaXast, mus, lambdas)
    meritproblem.M = problem0.M;
    meritproblem.cost = @(x) loneMeritFunction(problem0, x, rho);
    f0 = meritproblem.cost(xCur);
    df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
        problem0, mus, lambdas), deltaXast);
end