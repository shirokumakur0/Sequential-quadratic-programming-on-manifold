%%
%-------------------------- Fixed rank (Prototype) ------------------------------------
% NOT available for exactpenaltyViaMinimax (the solver doesn't support fixedrankembeddedfactory)

close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

row_dim = [2, 3, 4, 5, 7];
col_ratio = [1.5, 2, 2.5];
fixedrank_ratio = [0.25, 0.5, 0.75];
lc_dim_ratio = [0.1, 0.25, 0.5];

n_repeat = 1;   %Number of repeat on same set of data

for repeat = 1 : n_repeat
    for rdim = row_dim
        for cratio = col_ratio
            for rratio = fixedrank_ratio
                for dimratio = lc_dim_ratio
                    %_______Set up data______
                    cdim = ceil(cratio * rdim);
                    rank = ceil(rratio * rdim);
                    lc_dim = ceil((rdim + cdim) * dimratio);
                    
                    %______Set Object matrix_____
                    % Generate a random mxn matrix A of rank k
                    L = randn(rdim, rank);
                    R = randn(cdim, rank);
                    A = L*R';
                    % Generate a random mask for observed entries: P(i, j) = 1 if the entry
                    % (i, j) of A is observed, and 0 otherwise.
                    fraction = 4 * rank*(rdim+cdim-rank)/(rdim*cdim);
                    P = sparse(rand(rdim, cdim) <= fraction);
                    % Hence, we know the nonzero entries in PA:
                    PA = P.*A;
                    
                    %________Experiment_____
                    options.maxOuterIter = 1000;
                    options.maxtime = 3600;
                    options.minstepsize = 1e-8;
                    options.mineigval_correction = 1e-5;
                    options.verbosity = 1;

                    %________Setting________
                    setting.repeat = repeat;
                    setting.row_dim = rdim;
                    setting.col_dim = cdim;
                    setting.rank = rank;
                    setting.lc_dim = lc_dim;
                    setting.maxOuterIter = options.maxOuterIter;
                    setting.maxtime = options.maxtime;
                    setting.minstepsize = options.minstepsize;
                    setting.mineigval_correction = options.mineigval_correction;
                    setting.verbosity = options.verbosity;
                    setting.P = P;
                    setting.PA = PA;

                    %Only do mini-sum-max for low dimensional data
                    %if rdim == row_dim(1)
                    %    specifier.ind = ones(7,1);
                    %else
                    %    specifier.ind = [0, 1, 1, 1, 1, 1, 1];
                    %end

                    specifier.ind = [0,1,1,1,1,1,1];

                    result = clientconstraint_rank_constraints_nnlc_with_SQP(rdim, cdim, rank, lc_dim, P, PA, options, specifier, setting);
                    result = result(:);
                    param = [rdim; cdim; rank];
                    outputdata = [result; param]';

                    filename = sprintf('with_SQP_zz_RC_nnlc_RDim%dCDim%d.dat', rdim,cdim);
                    dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                end
            end
        end
    end
end