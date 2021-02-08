%%
%-------------------------- Fixed rank (Prototype) ------------------------------------
% NOT available for exactpenaltyViaMinimax (the solver doesn't support fixedrankembeddedfactory)

close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

tolKKTrespowerset = [16]; % 1e-* tolerance 
n_repeat = 1;   %Number of repeat on same set of data
minrank = 2;
for repeat = 1 : n_repeat
    for rdim = [4,5,6,7,10,15,20,25]  % to be [4,5,6,7,10,15,20,25]
        cdim = rdim * 2;
        maxrank = 3;
        % maxrank = max(ceil(0.5 * rdim), minrank);
        for tolKKTres = tolKKTrespowerset
            for rank = minrank:maxrank
                %_______Set up data______

                %______Set Object matrix_____
                % Generate a random mxn matrix A of rank k
                L = rand(rdim, rank);  % changed from randn to rand
                R = rand(cdim, rank);  % changed from randn to rand
                A =  L*R';
                % Generate a random mask for observed entries: P(i, j) = 1 if the entry
                % (i, j) of A is observed, and 0 otherwise.
                fraction = rank*(rdim+cdim-rank)/(rdim*cdim);  % 4* -> 1*, changed by MO
                P = sparse(rand(rdim, cdim) <= fraction);
                % Hence, we know the nonzero entries in PA:
                PA = P.*A;

                %________Experiment_____
                options.maxOuterIter = 100000;  % for Riemannian methods
                options.maxiter = options.maxOuterIter;  % for RSQP
                options.maxtime = 600;  % to be 600
                options.verbosity = 1;  % to be 1
                options.tolKKTres = 10^(-tolKKTres);     
                options.startingtolgradnorm = max(1e-3, 10^(-tolKKTres + 3));
                options.endingtolgradnorm = 10^(-tolKKTres);
                options.outerverbosity = options.verbosity;
                options.mineigval_correction = 1e-5;

                %________for initial point_____
                setting.initialpoint =  'feasible';
                %setting.initialpoint =  'random';

                %________Setting________
                setting.repeat = repeat;
                setting.row_dim = rdim;
                setting.col_dim = cdim;
                setting.rank = rank;
                setting.mineigval_correction = options.mineigval_correction;
                setting.tolKKTres =  tolKKTres;
                setting.maxOuterIter = options.maxOuterIter;
                setting.maxtime = options.maxtime;
                setting.verbosity = options.verbosity;
                setting.P = P;
                setting.PA = PA;

                specifier.ind = [1,1,1,1];  % to be [1,1,1,1]

                result = clientconstraint_rank_constraints_nn_with_SQP(rdim, cdim, rank, P, PA, options, specifier, setting);
                result = result(:);
                param = [rdim; cdim; rank; tolKKTres];
                outputdata = [result; param]';

                filename = sprintf('with_SQP_zz_RC_nn_Dim%dCDim%dTol%d.dat', rdim, cdim, tolKKTres);
                dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
            end
        end
    end
end