clear;
% Test 1, sphere case, where we use checkGradent function. Since this procedure
% generate a figure, if you don't want to make additional tabs,(for example, 
% when you run all tests at once?), please comment out this test.
problem = exampleProblemSphereHessian();
basis = makeCanonicalBasis(problem);
mus = ones(2,1);
lambdas = ones(1,1);
xCur = problem.M.rand();
subproblem = makeSubProbInfo(problem, xCur, mus, lambdas);
% Make the grad and hess of the Lagrangian at the current point
H = innerHess2matrix(subproblem.quadcost, basis);
%subproblem.hess = @(x, d) H*d;
hessMatLag =  hessMatLagrangian(xCur, subproblem, basis);
subproblem.hess = @(x, d) hessMatLag*d;
checkhessian(subproblem)
