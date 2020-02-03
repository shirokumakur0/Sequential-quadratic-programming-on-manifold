function newproblem = makeLagrangianforTest(problem, mus, lambdas)
% This is the function where we make a new subproblem whose obljective is
% the Lagrange function at the arguments.
% Note that this function is for tests, not for the original program (this
% is the reason why this function is at testdrivendevelopment) because in
% sqp.m, we don't need to make subproblem perfectly in actual since we
% finally convert the problem to what is suited for quadprog, which is very
% different from manopt structure.

% If you would make a subproblem as a struct (when you applied some Riemannian 
% method for subproblems?), please consider using tangentspacefactory.m,
% which is a well-established function on manopt.

costLag = @(X) costLagrangian(X, problem, mus, lambdas); % value
gradLag = @(X) gradLagrangian(X, problem, mus, lambdas); % in the tangent space
hessLag = @(X, d) hessLagrangian(X, d, problem, mus, lambdas); % in the tangent space
newproblem.cost = costLag;
newproblem.grad = gradLag;
newproblem.hess = hessLag;
newproblem.M = problem.M;

end