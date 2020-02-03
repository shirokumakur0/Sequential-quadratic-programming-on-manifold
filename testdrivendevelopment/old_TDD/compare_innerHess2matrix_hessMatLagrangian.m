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

primcostLag = @(X) costLagrangian(X, problem, mus, lambdas); % value
primgradLag = @(X) gradLagrangian(X, problem, mus, lambdas); % in the tangent space
primhessLag = @(X, d) hessLagrangian(X, d, problem, mus, lambdas); % in the tangent space
primproblem.cost = primcostLag;
primproblem.grad = primgradLag;
primproblem.hess = primhessLag;
primproblem.M = problem.M;
hessMatLag =  hessMatLagrangian(xCur, primproblem, basis);


% Test 2, stiefel case
problem = exampleProblemStiefelHessian();
mus = ones(2,1);
lambdas = ones(1,1);


% Test 3, oblique case
problem = exampleProblemObliqueHessian();
mus = ones(6,1);
lambdas = ones(1,1);
;