clear;
% This is a test for hessLagrangian.m, where we make a directional 
% hessian of the Lagrangian based on (cost/constraints) functions 
% with multipliers. The returned value 'hess' is a directional hessian of 
% the Lagrangian at x along d, which should be a vector.

% Test 1, sphere case
% manifold = spherefactory(2);
% problem.M = manifold;
% A = [1,2;2,1];
% problem.cost = @(x) x' * A * x;
% problem.egrad = @(x) 2*A*x ;
% problem.ehess = @(x,d) = 2*A*u;
% 
% I1 = [-1; 0];
% I2 = [0; 1];
% problem.ineq_constraint_cost = cell(2,1);
% problem.ineq_constraint_grad = cell(2,1);
% problem.ineq_constraint_cost{1} = @(x) I1' * x;
% problem.ineq_constraint_grad{1} = @(x) I1;
% problem.ineq_constraint_cost{2} = @(x) I2' * x;
% problem.ineq_constraint_grad{2} = @(x) I2;
% 
% E1 = [1;1];
% problem.eq_constraint_cost = cell(1,1);
% problem.eq_constraint_grad = cell(1,1);
% problem.eq_constraint_cost{1} = @(x) E1' * x;
% problem.eq_constraint_grad{1} = @(x) E1;
% 
% problem.condet = constraintsdetail(problem);
% mus = ones(2,1);
% lambdas = ones(1,1);
% % case 1
% x = [1;1];
% grad = gradLagrangian(x, problem, mus, lambdas);
% assert(isequal(grad,[-4;-1]));
% % case 2
% x = [0;1];
% grad = gradLagrangian(x, problem, mus, lambdas);
% assert(isequal(grad,[1;0]));
% % case 3
% x = [1;0];
% grad = gradLagrangian(x, problem, mus, lambdas);
% assert(isequal(grad,[0;4]));
