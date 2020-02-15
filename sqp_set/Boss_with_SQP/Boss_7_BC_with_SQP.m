%%
%--------------------------Balanced Cut------------------------------------
close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

%dim_set = [10, 25, 30, 50, 60]; %dimension of the Adjacency Matrix
dim_set = [5, 50, 75, 90, 100, 125, 150, 200];
density_set = [0.005, 0.01, 0.02, 0.04, 0.08]; %density of the Adjacency Matrix 


n_repeat = 1;   %Number of repeat on same set of data
rank = 2;     %Graph Bisection
seed_size = 5; %fixed seed size for BA
prob_ER = 0.5; %probability of connecting an edge in ER graph.

for repeat = 1 : n_repeat
    
    for dim = dim_set

        for density = density_set
        
            %_______Set up data______
            mlink = ceil(density * dim);
            L = powerlawgraph(seed_size, prob_ER, dim, mlink);
                        
            %________Experiment_____
            options.maxOuterIter = 5000;
            options.maxtime = 3600;
            options.minstepsize = 1e-4;
            options.verbosity = 1;

            %________Setting________
            setting.repeat = repeat;
            setting.dim = dim;
            setting.density = density;
            setting.mlink = mlink;
            setting.maxOuterIter = options.maxOuterIter;
            setting.maxtime = options.maxtime;
            setting.minstepsize = options.minstepsize;
            setting.verbosity = options.verbosity;
            %setting.trimhess_perturbation = options.trimhess_perturbation;
            setting.L = L;
            
            
            %Only do mini-sum-max for low dimensional data
            if dim == dim_set(1)
                specifier.ind = ones(7,1);
            else
                specifier.ind = [0, 1, 1, 1, 1, 1, 1];
            end
            
            result = clientconstraint_oblique_balancedcut_with_SQP(L, rank, options, specifier, setting);
            result = result(:);
            param = [dim; density; repeat];
            outputdata = [result; param]';
            
            filename = sprintf('with_SQP_zz_BC_Dim%d.dat', dim);
            dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
        end
        
    end
    
end
