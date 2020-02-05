    clear;
    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt.']);
    end
    
    % Generate the problem data.
    n = 3;
    A = zeros(n);
    A(n,n) = 10;
    A(1,1) = 1;
    A(2,2) = 100;
    
    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Define the problem cost function and its gradient.
    problem.cost  = @(x) x'*(A*x);
    problem.egrad = @(x) 2*A*x;
    problem.ehess = @(x, xdot) 2*A*xdot;
    
    I1 = [-1; 0; 0];
    I2 = [0; -1; 0];
    I3 = [0; 0; -1];

    problem.ineq_constraint_cost = cell(3,1);
    problem.ineq_constraint_grad = cell(3,1);
    problem.ineq_constraint_hess = cell(3,1);
    problem.ineq_constraint_cost{1} = @(x) I1' * x;
    problem.ineq_constraint_grad{1} = @(x) I1;
    problem.ineq_constraint_hess{1} = @(x, u) 0;
    problem.ineq_constraint_cost{2} = @(x) I2' * x;
    problem.ineq_constraint_grad{2} = @(x) I2;
    problem.ineq_constraint_hess{2} = @(x, u) 0;
    problem.ineq_constraint_cost{3} = @(x) I3' * x;
    problem.ineq_constraint_grad{3} = @(x) I3;
    problem.ineq_constraint_hess{3} = @(x, u) 0;
    
    % Numerically check gradient and Hessian consistency.
    %     figure;
    %     checkgradient(problem);
    %     figure;
    %     checkhessian(problem);

    % Solve.
    % if you use ALM    
%          options.maxOuterIter = 10000000;
%          options.maxtime = 3600;
%          options.minstepsize = 1e-10;
     x0 = problem.M.rand();
%          [xfinal, info] = almbddmultiplier(problem, x0, options);
    options.trimhessian = "mineigval_manopt";
    % if you use sqp
     [x, xcost, info] = sqp(problem, x0, options); % #ok<ASGLU>    
    
    % Display some statistics.
    figure;
    semilogy([info.iter], [info.cost], '.-');
    xlabel('Iteration #');
    ylabel('cost');
    title('Convergence of the sqp algorithm on the nonnegative sphere');
