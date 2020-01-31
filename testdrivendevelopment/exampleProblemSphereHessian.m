function problem = exampleProblemSphereHessian()
% This is the example for tests which we will use to verify SQP programs.
% The followig manifold is easy and hence we can understand the properties
% of this manifold (and outputs from functions when input this one). Here
% we consider a sphere case with 2 ineqalities and 1 equality constraints.

% The difference between exampleProblemSphere and this file is the
% Hessian functions of the objective and constraints functions. In this
% time, we changed the objective funcions and constraints to the quadratic
% ones to make these of Hessians untrivial.

manifold = spherefactory(2);
problem.M = manifold;

A = rand(2,2);
problem.cost = @(X) X' * A * X;
problem.egrad = @(X) A*X +A'*X;
problem.ehess = @(X, U) A*U +A'*U;

I1 = rand(2,2);
I2 = rand(2,2);

problem.ineq_constraint_cost = cell(2,1);
problem.ineq_constraint_grad = cell(2,1);
problem.ineq_constraint_hess = cell(2,1);
problem.ineq_constraint_cost{1} = @(X) X' * I1 * X;
problem.ineq_constraint_grad{1} = @(X) I1*X +I1'*X;
problem.ineq_constraint_hess{1} = @(X, U) I1*U +I1'*U;
problem.ineq_constraint_cost{2} = @(X) X' * I2 * X;
problem.ineq_constraint_grad{2} = @(X) I2*X +I2'*X;
problem.ineq_constraint_hess{2} = @(X, U) I2*U +I2'*U;

E1 = rand(2,2);
problem.eq_constraint_cost = cell(1,1);
problem.eq_constraint_grad = cell(1,1);
problem.eq_constraint_hess = cell(1,1);
problem.eq_constraint_cost{1} = @(X) X' * E1 * X;
problem.eq_constraint_grad{1} = @(X) E1*X+E1'*X;
problem.eq_constraint_hess{1} = @(X, U) E1*U+E1'*U;
problem.condet = constraintsdetail(problem);

end