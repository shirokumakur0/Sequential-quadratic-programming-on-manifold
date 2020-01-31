function problem = exampleProblemOblique()

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

end