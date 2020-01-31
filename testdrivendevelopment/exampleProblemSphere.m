function problem = exampleProblemSphere()
% This is the example for tests which we will use to verify SQP programs.
% The followig manifold is easy and hence we can understand the properties
% of this manifold (and outputs from functions when input this one). Here
% we consider a sphere case with 2 ineqalities and 1 equality constraints.

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

end