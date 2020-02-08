function localdefaults = setLocalDefaults(problem0)

    % Local defaults for the program
    localdefaults.maxiter = 1000;
    localdefaults.maxtime = 3600;
    localdefaults.minstepsize = 1e-6;
    localdefaults.tolgradnorm = 1e-7;
    localdefaults.storedepth = 3;
    localdefaults.tau = 0.8;  % TODO: should find an appropriate value as long as tau > 0
    localdefaults.rho = 1;  % TODO: should find an appropriate value as long as rho > 0
    localdefaults.beta = 0.5;  % TODO: should find an appropriate value as long as 1 > beta > 0
    localdefaults.gamma = 0.5; % TODO: should find an appropriate value as long as 1 > gamma > 0  
    localdefaults.mus = ones(problem0.condet.n_ineq_constraint_cost, 1);
    localdefaults.lambdas = ones(problem0.condet.n_eq_constraint_cost, 1);    
    localdefaults.ls_max_steps  = 30;
    localdefaults.regularhesseigval = 1e-3;
    localdefaults.trimhessian = 'mineigval_manopt';
end
