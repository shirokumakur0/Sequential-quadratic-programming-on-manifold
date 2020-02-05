function qpinfo = trimHessianMatrix(qpinfo, options)
    if strcmp(options.trimhessian, "eye")
        qpinfo.H = eye(qpinfo.n);
    else % elseif options.trimhessian == 'mineigval'
        eigval = eig(qpinfo.H);
        mineigval = min(eigval);
        fprintf("min eigen val from trimHessian: %f", mineigval);
        if mineigval < 0
            regularize_coeff = max(options.regularhesseigval, abs(mineigval));
            qpinfo.H = qpinfo.H + regularize_coeff * eye(qpinfo.n);
        end
    end
end