%%
%--------------------------Balanced Cut------------------------------------
close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

% (original)
%dim_set = [10, 25, 30, 50, 60]; %dimension of the Adjacency Matrix
%dim_set = [10, 50, 100, 200, 500, 1000];
%density_set = [0.005, 0.01, 0.02, 0.04, 0.08]; %density of the Adjacency Matrix 

dim_set = [10, 50, 75, 100, 250];
density_set = [0.005 ,0.01, 0.02, 0.03, 0.04];
tolKKTrespowerset = [2, 4, 6, 8, 10]; % 1e-* tolerance

n_repeat = 1;   %Number of repeat on same set of data
rank = 2;     %Graph Bisection
seed_size = 5; %fixed seed size for BA
prob_ER = 0.5; %probability of connecting an edge in ER graph.

for repeat = 1 : n_repeat
    
    for dim = dim_set

        for density = density_set
            
            for tolKKTres = tolKKTrespowerset

                %_______Set up data______
                mlink = ceil(density * dim);
                L = powerlawgraph(seed_size, prob_ER, dim, mlink);

                %________Experiment_____
                options.maxOuterIter = 10000;  % for Riemannian methods & fmincon
                options.maxiter = options.maxOuterIter;  % for RSQP
                options.maxtime = 180;
                options.tolKKTres = 10^(-tolKKTres);
                %options.beta = 0.5; % for SQP
                options.verbosity = 1;
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

                specifier.ind = [1,1,1,1,1,1];
                
                result = clientconstraint_modified_oblique_balancedcut_with_SQP(L, rank, options, specifier, setting);
                result = result(:);
                param = [dim; density; repeat; tolKKTres];
                outputdata = [result; param]';
                
                % PP according to dimension and residual
                filename = sprintf('with_SQP_zz_BC_Dim%dTol%d.dat', dim,tolKKTres);
                dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                
                % PP according to dimension
                %filename = sprintf('with_SQP_zz_BC_Dim%d.dat', dim);
                %dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                % PP according to dimension
                %filename = sprintf('with_SQP_zz_BC_Tol%d.dat', tolKKTres);
                %dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
            end
        end
    end
    
end
