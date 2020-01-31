clear;
% This is a test for costLagrangian.m where we make a cost of the Lagrangian
% based on values (cost/constraints) functions values with multipliers.
% The returned value 'val' is the scalar of the Lagrangian at x. 

% Test 1, sphere case
problem = exampleProblemSphere();
mus = ones(2,1);
lambdas = ones(1,1);
% case 1
x = [1;1];
val = costLagrangian(x, problem, mus, lambdas);
assert(val == 5);
% case 2
x = [2;1];
val = costLagrangian(x, problem, mus, lambdas);
assert(val == 6);
% case 3
x = [0;1];
val = costLagrangian(x, problem, mus, lambdas);
assert(val == 4);
% case 4
x = [1;0];
val = costLagrangian(x, problem, mus, lambdas);
assert(val == 1)

% Test 2, Stiefel
problem = exampleProblemStiefel();
mus = ones(6,1);
lambdas = ones(1,1);
% case 1
X = [1,0;1,0;1,0];
val = costLagrangian(X, problem, mus, lambdas);
assert(val == 12);
% case 2
X = [0,1;0,1;0,1];
val = costLagrangian(X, problem, mus, lambdas);
assert(val == 0);

% Test 3, oblique
problem = exampleProblemOblique();
mus = ones(6,1);
lambdas = ones(1,1);
% case 1
X = [1,1,1;0,0,0];
val = costLagrangian(X, problem, mus, lambdas);
assert(val ==12);
% case 2
X = [0,0,0;1,1,1];
val = costLagrangian(X, problem, mus, lambdas);
assert(val == 21);

fprintf('All tests have been accepted! [test_costLagrangian]\n')
