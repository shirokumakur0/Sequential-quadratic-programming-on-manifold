%%
%-------------------------- Fixed rank with nonnegativity and equality ----
%-------------------------- Measurement of Accuracy -----------------------
close all; clear all; clc;
specifier.matlabversion = 0; % 0 if older than 2015 1 otherwise

tolKKTrespowerset = [16]; % log10 scale, i.e., 1e-* tolerance 
n_repeat = 1;   % Number of repeat on same set of data
minrank = 2;
for repeat = 1 : n_repeat
    for rdim = [4,5,6,7]  % [4,5,6,7]
        cdim = rdim * 2;
        maxrank = minrank;  % or, e.g., maxrank = 3;
        for tolKKTres = tolKKTrespowerset
            for rankval = minrank:maxrank
                for eqconstratio = [0.5]  % [0.5]
                    for maskratio = [0.5]  % [0.5]
                        %_______Set up data______

                        %______Set Object matrix_____
                        % Generate a random mxn matrix A of rank k
                        while true
                            L = rand(rdim, rankval);
                            R = rand(cdim, rankval);
                            A =  L*R';
                            rankA = rank(A);
                            if rankA == rankval
                                break
                            end
                        end
                        % Generate a random mask for observed entries: P(i, j) = 1 if the entry
                        % (i, j) of A is observed, and 0 otherwise.
                        initP = zeros(rdim, cdim);
                        observed_indices = randperm(rdim*cdim, ceil(maskratio * rdim*cdim));
                        for obs = observed_indices
                            initP(obs) = 1;
                        end
                        
                        %%% For nonnegativity inequality and equality
                        nonzero_num = nnz(initP);
                        eqnum = ceil(nonzero_num * eqconstratio);
                        nonzeroidcs = find(initP);
                        s = RandStream('mlfg6331_64'); 
                        eqindices = sort(randsample(s,nonzeroidcs,eqnum));
                        for eqindex = eqindices'
                            initP(eqindex) = 0;
                        end
                        P = initP;
                        
                        %________Experiment_____
                        options.maxOuterIter = 100000;  % for Riemannian methods
                        options.maxiter = options.maxOuterIter;  % for RSQO (RSQP)
                        options.maxtime = 600;  % 600
                        options.verbosity = 1;  % 1
                        options.tolKKTres = 10^(-tolKKTres);     
                        options.startingtolgradnorm = max(1e-3, 10^(-tolKKTres + 3));
                        options.endingtolgradnorm = 10^(-tolKKTres);
                        options.outerverbosity = options.verbosity;
                        options.mineigval_correction = 1e-5;  % 1e-5

                        %________for initial point_____
                        setting.initialpoint =  'feasible_region';  % other cndidates: 'eye', 'random'
                        
                        %________Setting________
                        setting.repeat = repeat;
                        setting.row_dim = rdim;
                        setting.col_dim = cdim;
                        setting.rank = rankval;
                        setting.mineigval_correction = options.mineigval_correction;
                        setting.tolKKTres =  tolKKTres;
                        setting.maxOuterIter = options.maxOuterIter;
                        setting.maxtime = options.maxtime;
                        setting.verbosity = options.verbosity;
                        setting.eqconstratio = eqconstratio;
                        setting.maskratio = maskratio;
                        setting.filepath = sprintf('nrep%dRowdim%dColdim%dRank%dTol%dEqratio%0.1eMaskratio%0.1e',...
                            setting.repeat, setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres,setting.eqconstratio, setting.maskratio);
                        setting.P = P;
                        setting.A = A;

                        specifier.ind = [1,1,1,1];  % [1,1,1,1]

                        result = clientconstraint_fixedrank_equality_nonneg_lowrankcompletion(rdim, cdim, rankval, P, A, eqindices, options, specifier, setting);
                        result = result(:);
                        param = [rdim; cdim; rankval; tolKKTres];
                        outputdata = [result; param]';
                        
                        % Performance profile
                        filename = sprintf('with_SQP_zz_RC_nnlc_%s.dat', setting.filepath);
                        dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                    end
                end
            end
        end
    end
end