    clear;
    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt.']);
    end
    
    % Generate the problem data.
    n = 2;
    A = randn(n);
    A = .5*(A+A');
    
    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Define the problem cost function and its gradient.
    problem.cost  = @(x) -x'*(A*x);
    problem.egrad = @(x) -2*A*x;
    problem.ehess = @(x, xdot) -2*A*xdot;
    
    I1 = [-1; 0];
    I2 = [0; 1];
    problem.ineq_constraint_cost = cell(2,1);
    problem.ineq_constraint_grad = cell(2,1);
    problem.ineq_constraint_hess = cell(2,1);
    problem.ineq_constraint_cost{1} = @(x) I1' * x;
    problem.ineq_constraint_grad{1} = @(x) I1;
    problem.ineq_constraint_hess{1} = @(x, u) 0;
    problem.ineq_constraint_cost{2} = @(x) I2' * x;
    problem.ineq_constraint_grad{2} = @(x) I2;
    problem.ineq_constraint_hess{2} = @(x, u) 0;

    E1 = [1;1];
    problem.eq_constraint_cost = cell(1,1);
    problem.eq_constraint_grad = cell(1,1);
    problem.eq_constraint_hess = cell(1,1);
    problem.eq_constraint_cost{1} = @(x) E1' * x;
    problem.eq_constraint_grad{1} = @(x) E1;
    problem.eq_constraint_hess{1} = @(x, u) 0;

    % problem.condet = constraintsdetail(problem);
    
    % Numerically check gradient and Hessian consistency.
    % figure;
    % checkgradient(problem);
    % figure;
    % checkhessian(problem);
 
    % Solve.
    %[x, xcost, info] = steepestdescent(problem);
    [x, xcost] = sqp(problem, [], []); % #ok<ASGLU>
    
    % Display some statistics.
    % figure;
    % semilogy([info.iter], [info.gradnorm], '.-');
    % xlabel('Iteration #');
    % ylabel('Gradient norm');
    % title('Convergence of the steepest-descent algorithm on the sphere');
