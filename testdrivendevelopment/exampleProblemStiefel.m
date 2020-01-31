function problem = exampleProblemStiefel()
% This is the example for tests which we will use to verify SQP programs.
% The followig manifold is easy and hence we can understand the properties
% of this manifold (and outputs from functions when input this one). Here
% we consider a stiefel case with 2 ineqalities and 1 equality constraints.

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

end