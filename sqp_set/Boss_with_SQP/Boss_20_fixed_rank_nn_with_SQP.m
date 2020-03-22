%%
%-------------------------- Fixed rank (Prototype) ------------------------------------
% NOT available for exactpenaltyViaMinimax (the solver doesn't support fixedrankembeddedfactory)

close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

row_dim =[2,3,4];
%col_ratio = [1.5, 1.75 ,2];
% fixedrank_ratio = [0.1, 0.25, 0.5];
%tolKKTrespowerset = [1, 2, 3, 4, 5, 6]; % 1e-* tolerance

n_repeat = 8;   %Number of repeat on same set of data

for repeat = 1 : n_repeat
    for rdim = row_dim
        fixedrank_repeat = ceil(rdim / 2);
        tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance
        %if rdim == 7
        %    tolKKTrespowerset = [2, 3, 4, 5]; % 1e-* tolerance
        %elseif rdim >= 6
        %    tolKKTrespowerset = [2, 3, 4, 5, 6]; % 1e-* tolerance
        %else
        %    tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance
        %end
      
        if (rdim == 2) || (rdim == 3)
            col_ratio = [1.5, 2, 3];
        else
            col_ratio = [1.5, 1.75, 2, 3];
        end
        
        for cratio = col_ratio
            % for rratio = fixedrank_ratio
            for rank = 1:fixedrank_repeat
                for tolKKTres = tolKKTrespowerset
                    %_______Set up data______
                    cdim = ceil(cratio * rdim);
                    %rank = ceil(rratio * rdim);
 
                    %______Set Object matrix_____
                    % Generate a random mxn matrix A of rank k
                    L = randn(rdim, rank);
                    R = randn(cdim, rank);
                    A = L*R';
                    % Generate a random mask for observed entries: P(i, j) = 1 if the entry
                    % (i, j) of A is observed, and 0 otherwise.
                    fraction = rank*(rdim+cdim-rank)/(rdim*cdim);  % 4* -> 1*, changed by MO
                    P = sparse(rand(rdim, cdim) <= fraction);
                    % Hence, we know the nonzero entries in PA:
                    PA = P.*A;
                    
                    %________for initial point_____
                    setting.initialpoint =  "feasible";
                    %setting.initialpoint =  "random";
                        
                    
                    %________Experiment_____
                    options.maxOuterIter = 10000;  % for Riemannian methods
                    options.maxiter = options.maxOuterIter;  % for RSQP
                    options.maxtime = 180;
                    options.verbosity = 1;
                    options.tolKKTres = 10^(-tolKKTres);                    
                    options.outerverbosity = options.verbosity;
                    options.mineigval_correction = 1e-5;

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

                    specifier.ind = [1,1,1,1];

                    result = clientconstraint_rank_constraints_nn_with_SQP(rdim, cdim, rank, P, PA, options, specifier, setting);
                    result = result(:);
                    param = [rdim; cdim; rank; tolKKTres];
                    outputdata = [result; param]';
                    
                    filename = sprintf('with_SQP_zz_RC_nn_RDim%dCDim%dTol%d.dat', rdim, cdim, tolKKTres);
                    dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                    
                    % PP according to dimension
                    %filename = sprintf('with_SQP_zz_RC_nn_RDim%dCDim%d.dat', rdim,cdim);
                    %dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                    % PP according to tolKKTres
                    %filename = sprintf('with_SQP_zz_RC_nn_Tol%d.dat', tolKKTres);
                    %dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                end
            end
        end
    end
end