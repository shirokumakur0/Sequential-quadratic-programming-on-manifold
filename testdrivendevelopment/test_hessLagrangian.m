clear;
% This is a test for hessLagrangian.m, where we make a directional 
% hessian of the Lagrangian based on (cost/constraints) functions 
% with multipliers. The returned value 'hess' is a directional hessian of 
% the Lagrangian at x along d, which should be a vector.

% Test 1, sphere case, where we use checkGradent function. Since this procedure
% generate a figure, if you don't want to make additional tabs,(for example, 
% when you run all tests at once?), please comment out this test.
problem = exampleProblemSphereHessian();
mus = ones(2,1);
lambdas = ones(1,1);
newproblem = makeLagrangianforTest(problem,mus,lambdas);
checkhessian(newproblem);

% Test 2, stiefel case
problem = exampleProblemStiefelHessian();
mus = ones(2,1);
lambdas = ones(1,1);
newproblem = makeLagrangianforTest(problem,mus,lambdas);
checkhessian(newproblem);

% Test 3, oblique case
problem = exampleProblemObliqueHessian();
mus = ones(6,1);
lambdas = ones(1,1);
newproblem = makeLagrangianforTest(problem,mus,lambdas);
checkhessian(newproblem);

fprintf('All tests have done! Check the generated figures! [test_hessLagrangian]\n')