clear;
% This is a test for costLagrangian.m where we make a cost of the Lagrangian
% based on values (cost/constraints) functions values with multipliers.
% The returned value 'val' is the scalar of the Lagrangian at x. 

% Test 1, sphere case
manifold = spherefactory(2);
problem.M = manifold;
A = [1;2];
problem.cost = @(x) A' * x;
problem.egrad = @(x) A;

I1 = [-1; 0];
I2 = [0; 1];
problem.ineq_constraint_cost = cell(2,1);
problem.ineq_constraint_grad = cell(2,1);
problem.ineq_constraint_cost{1} = @(x) I1' * x;
problem.ineq_constraint_grad{1} = @(x) I1;
problem.ineq_constraint_cost{2} = @(x) I2' * x;
problem.ineq_constraint_grad{2} = @(x) I2;

E1 = [1;1];
problem.eq_constraint_cost = cell(1,1);
problem.eq_constraint_grad = cell(1,1);
problem.eq_constraint_cost{1} = @(x) E1' * x;
problem.eq_constraint_grad{1} = @(x) E1;

problem.condet = constraintsdetail(problem);
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
manifold = stiefelfactory(3,2);
problem.M = manifold;
A = [1,-1;2,-2;3,-3];
problem.cost = @(X) trace(A' * X);
problem.egrad = @(X) A;

I1 = zeros(3,2);
I1(1) = 1;
I2 = zeros(3,2);
I2(2) = 1;
I3 = zeros(3,2);
I3(3) = 1;
I4 = zeros(3,2);
I4(4) = 1;
I5 = zeros(3,2);
I5(5) = 1;
I6 = zeros(3,2);
I6(6) = 1;
problem.ineq_constraint_cost = cell(6,1);
problem.ineq_constraint_grad = cell(6,1);
problem.ineq_constraint_cost{1} = @(X) trace(I1' * X);
problem.ineq_constraint_grad{1} = @(X) I1;
problem.ineq_constraint_cost{2} = @(X) trace(I2' * X);
problem.ineq_constraint_grad{2} = @(X) I2;
problem.ineq_constraint_cost{3} = @(X) trace(I3' * X);
problem.ineq_constraint_grad{3} = @(X) I3;
problem.ineq_constraint_cost{4} = @(X) trace(I4' * X);
problem.ineq_constraint_grad{4} = @(X) I4;
problem.ineq_constraint_cost{5} = @(X) trace(I5' * X);
problem.ineq_constraint_grad{5} = @(X) I5;
problem.ineq_constraint_cost{6} = @(X) trace(I6' * X);
problem.ineq_constraint_grad{6} = @(X) I6;

E1 = [1,1;1,1;1,1];
problem.eq_constraint_cost = cell(1,1);
problem.eq_constraint_grad = cell(1,1);
problem.eq_constraint_cost{1} = @(X) trace(E1' * X);
problem.eq_constraint_grad{1} = @(X) E1;

problem.condet = constraintsdetail(problem);
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
manifold = obliquefactory(2,3);
problem.M = manifold;
A = [1,2,3;4,5,6];
problem.cost = @(X) trace(A' * X);
problem.egrad = @(X) A;

I1 = zeros(2,3);
I1(1) = 1;
I2 = zeros(2,3);
I2(2) = 1;
I3 = zeros(2,3);
I3(3) = 1;
I4 = zeros(2,3);
I4(4) = 1;
I5 = zeros(2,3);
I5(5) = 1;
I6 = zeros(2,3);
I6(6) = 1;
problem.ineq_constraint_cost = cell(6,1);
problem.ineq_constraint_grad = cell(6,1);
problem.ineq_constraint_cost{1} = @(X) trace(I1' * X);
problem.ineq_constraint_grad{1} = @(X) I1;
problem.ineq_constraint_cost{2} = @(X) trace(I2' * X);
problem.ineq_constraint_grad{2} = @(X) I2;
problem.ineq_constraint_cost{3} = @(X) trace(I3' * X);
problem.ineq_constraint_grad{3} = @(X) I3;
problem.ineq_constraint_cost{4} = @(X) trace(I4' * X);
problem.ineq_constraint_grad{4} = @(X) I4;
problem.ineq_constraint_cost{5} = @(X) trace(I5' * X);
problem.ineq_constraint_grad{5} = @(X) I5;
problem.ineq_constraint_cost{6} = @(X) trace(I6' * X);
problem.ineq_constraint_grad{6} = @(X) I6;

E1 = [1,1,1;1,1,1];
problem.eq_constraint_cost = cell(1,1);
problem.eq_constraint_grad = cell(1,1);
problem.eq_constraint_cost{1} = @(X) trace(E1' * X);
problem.eq_constraint_grad{1} = @(X) E1;

problem.condet = constraintsdetail(problem);
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