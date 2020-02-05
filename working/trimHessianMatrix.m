function qpinfo = trimHessianMatrix(qpinfo, options)
    if strcmp(options.trimhessian, "eye")
        qpinfo.H = eye(qpinfo.n);
    elseif strcmp(options.trimhessian, 'mineigval_matlab')
        eigval = eig(qpinfo.H);
        mineigval = min(eigval);
        if mineigval < 0
            regularize_coeff = max(options.regularhesseigval, abs(mineigval));
            qpinfo.H = qpinfo.H + regularize_coeff * eye(qpinfo.n);
        end
    elseif strcmp(options.trimhessian, 'mineigval_manopt')
        % the difference between mineiegval_manopt and mineigval_matlab is
        % a solver for calculating the minimum eigen value. _matlab use the
        % eigenvalue decomposition function on matlab, whereas, _manopt use
        % the hessianextreme function which formulzize and solve
        % an optimization problem to get a minimum eigenvalue and 
        % the corresponding eigenvector.
        % However, here we just get a minimum eigenval which must be
        % calculated in makeQPInfo.m.
        mineigval = qpinfo.mineigval_manopt;
        if mineigval < 0
            regularize_coeff = max(options.regularhesseigval, abs(mineigval));
            qpinfo.H = qpinfo.H + regularize_coeff * eye(qpinfo.n);
        end
    end
end