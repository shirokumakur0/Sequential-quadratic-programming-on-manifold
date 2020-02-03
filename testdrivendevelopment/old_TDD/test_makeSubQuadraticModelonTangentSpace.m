clear;
% This is a test for makeSubQuadraticModelonTangentSpace.m, where we make a
% subproblem with the Lagrangian and other constraints from the original
% problem (and situatitons, that is, Lagrange multipliers).
% We assume that solvers for the ganerated subproblems are first-order
% methods, thus, we don't construct any hessian function for objective or
% constraints functions.

% Test 1, sphere case, where we use checkGradent function. Since this procedure
% generate a figure, if you don't want to make additional tabs,(for example, 
% when you run all tests at once?), please comment out this test.
problem = exampleProblemSphereHessian();
xCur = problem.M.rand();
mus = ones(2,1);
lambdas = ones(1,1);
subproblem = makeSubQuadraticModelonTangentSpace(problem, xCur, mus, lambdas);
figure;
checkgradient(subproblem);

xTan = subproblem.M.rand();
dTan = subproblem.M.randvec(xTan);

figure;
checkconstraints(subproblem,xTan,dTan);

% Test 2, stiefel case
problem = exampleProblemStiefelHessian();
xCur = problem.M.rand();
mus = ones(2,1);
lambdas = ones(1,1);
subproblem = makeSubQuadraticModelonTangentSpace(problem, xCur, mus, lambdas);
% figure;
% checkgradient(subproblem);

xTan = subproblem.M.rand();
dTan = subproblem.M.randvec(xTan);

figure;
checkconstraints(subproblem,xTan,dTan);

% Test 3, oblique case
problem = exampleProblemObliqueHessian();
mus = ones(6,1);
lambdas = ones(1,1);
xCur = problem.M.rand();
subproblem = makeSubQuadraticModelonTangentSpace(problem, xCur, mus, lambdas);
checkgradient(subproblem);
xTan = subproblem.M.rand();
dTan = subproblem.M.randvec(xTan);

figure;
checkconstraints(subproblem,xTan,dTan);

fprintf('All tests have done! Check the generated figures! [test_hessLagrangian]\n')