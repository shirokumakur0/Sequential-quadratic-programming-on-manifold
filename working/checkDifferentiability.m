function problem0 = checkDifferentiability(problem0)
% The funcion is imported from a part of an usual solver on manopt, as it is.
    if ~canGetCost(problem0)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    if ~canGetGradient(problem0) && ~canGetApproxGradient(problem0)
        % Note: we do not give a warning if an approximate gradient is
        % explicitly given in the problem description, as in that case the user
        % seems to be aware of the issue.
        warning('manopt:getGradient:approx', ...
           ['No gradient provided. Using an FD approximation instead (slow).\n' ...
            'It may be necessary to increase options.tolgradnorm.\n' ...
            'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
        problem0.approxgrad = approxgradientFD(problem0);
   end
end