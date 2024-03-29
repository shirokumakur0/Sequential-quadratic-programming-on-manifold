%%
%--------------------------Balanced Cut------------------------------------
close all; clear all; clc;
specifier.matlabversion = 0; % 0 if older than 2015 1 otherwise

dim_set = [50];
density_set = [0.01];
tolKKTrespowerset = [16]; % log10 scale, i.e., 1e-* tolerance

n_repeat = 1;   % Number of repeat on same set of data
rank = 2;     % Graph Bisection
seed_size = 5; % fixed seed size for BA
prob_ER = 0.5; % probability of connecting an edge in ER graph.

for repeat = 1 : n_repeat
    
    for dim = dim_set

        for density = density_set
            
            for tolKKTres = tolKKTrespowerset

                %_______Set up data______
                mlink = ceil(density * dim);
                L = powerlawgraph(seed_size, prob_ER, dim, mlink);

                %________Experiment_____
                options.maxOuterIter = 100000;  % for Riemannian methods & fmincon
                options.maxiter = options.maxOuterIter;  % for RSQO (RSQP)
                options.maxtime = 600;  % 600
                options.tolKKTres = 10^(-tolKKTres);
                options.startingtolgradnorm = max(1e-3,10^(-tolKKTres + 3));
                options.endingtolgradnorm = 10^(-tolKKTres);
                options.verbosity = 1;  % 1
                options.outerverbosity = options.verbosity;

                %________Setting________
                setting.repeat = repeat;
                setting.dim = dim;
                setting.density = density;
                setting.mlink = mlink;
                setting.tolKKTres = tolKKTres;
                setting.maxOuterIter = options.maxOuterIter;
                setting.maxtime = options.maxtime;
                setting.L = L;

                specifier.ind = [1,1,1,1,1];  % [1,1,1,1,1]
                
                result = clientconstraint_oblique_balancedcut(L, rank, options, specifier, setting);
                result = result(:);
                param = [dim; density; repeat; tolKKTres];
                outputdata = [result; param]';
                
                % Performance profile according to dimension and residual
                % filename = sprintf('with_SQP_zz_BC_Dim%dTol%d.dat', dim,tolKKTres);
                % dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                
            end
        end
    end
    
end
