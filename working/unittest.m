% Generate random problem data.
n = 1000;
A = randn(n);
A = .5*(A+A.');
 
% Create the problem structure.
manifold = spherefactory(n);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) -x'*(A*x);
problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean
 
% Numerically check gradient consistency (optional).
checkgradient(problem);
 
% Solve.
[x, xcost, info, options] = working_sqponmani(problem);
