clear;
% This is a test for gradLagrangian.m, where we make a gradient of the 
% Lagrangian based on values (cost/constraints) functions values
% with multipliers. The returned value 'grad' is a gradient of 
% the Lagrangian at x, which should be a vector.

% Test 1, sphere case
problem = exampleProblemSphere();
mus = ones(2,1);
lambdas = ones(1,1);
% case 1
x = [1;1];
grad = gradLagrangian(x, problem, mus, lambdas);
assert(isequal(grad,[-4;-1]));
% case 2
x = [0;1];
grad = gradLagrangian(x, problem, mus, lambdas);
assert(isequal(grad,[1;0]));
% case 3
x = [1;0];
grad = gradLagrangian(x, problem, mus, lambdas);
assert(isequal(grad,[0;4]));

% case 4, where we use checkGradent function. Since this procedure
% generate a figure, if you don't want to make additional tabs,(for example, 
% when you run all tests at once?), please comment out this test.

% newproblem = makeLagrangianforTest(problem,mus,lambdas);
% checkgradient(newproblem);

% Test 2, stiefel case
problem = exampleProblemStiefel();
mus = ones(6,1);
lambdas = ones(1,1);
% newproblem = makeLagrangianforTest(problem,mus,lambdas);
% checkgradient(newproblem);

% Test 3, oblique case
problem = exampleProblemOblique();
mus = ones(6,1);
lambdas = ones(1,1);
% newproblem = makeLagrangianforTest(problem,mus,lambdas);
% checkgradient(newproblem);

fprintf('All tests have been accepted! [test_gradLagrangian]\n')