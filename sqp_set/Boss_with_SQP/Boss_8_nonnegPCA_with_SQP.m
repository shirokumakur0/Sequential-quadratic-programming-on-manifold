%%
%---------------------------------PCA------------------------------------

close all; clc; clear all;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

dim_set = [10, 50, 100, 150, 250, 500, 1000];  % Dimension of "the Cov Matrix"
snrset = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]; % Signal Strength
deltaset = [0.1, 0.3, 0.7, 0.9];           % Sparsity 
rank = 1;                                  % Rank of BM Relaxation. 1 if we don't.
n_repeat = 1;                              % Number of repeat experiment

for repeat = 1: n_repeat
    
    for dim = dim_set
        
        for snr = snrset
            
            for delta = deltaset

                %_______Set up data______
                T = dim;
                samplesize = floor(delta*dim);
                S = randsample(dim, samplesize);
                v = zeros(dim,1);
                v(S) = 1/sqrt(samplesize);
                X = sqrt(snr) * v * (v.');
                Z = randn(dim)/sqrt(T);
                for ii = 1: dim
                    Z(ii,ii) = randn * 2/sqrt(T);
                end
                X = X+Z;
                

                %________Experiment_____
                options.maxOuterIter = 5000;
                options.maxtime = 3600;
                options.minstepsize = 1e-4;
                options.mineigval_correction = 0;
                options.mineigval_correction = 1e-10;
                options.verbosity = 1;
                
                %________Setting________
                setting.repeat = repeat;
                setting.dim = dim;
                setting.snr = snr;
                setting.delta = delta;
                setting.rank = rank;
                setting.maxOuterIter = options.maxOuterIter;
                setting.maxtime = options.maxtime;
                setting.minstepsize = options.minstepsize;
                setting.verbosity = options.verbosity;
                setting.Z = Z;
                
                %Only do mini-sum-max for low dimensional data
                if dim == dim_set(1)
                    specifier.ind = ones(7,1);
                else
                    specifier.ind = [0, 1, 1, 1, 1, 1, 1];
                end
                %specifier.ind = [0, 0, 0, 0, 0, 0, 1];

                result = clientconstraint_sphere_nonnegativePCA_with_SQP(X, rank, options, specifier, setting);
                result = result(:);
                param = [dim; snr; delta; repeat];
                outputdata = [result; param]';
                
                filename = sprintf('with_SQP_zz_NNPCA_Dim%d.dat', dim);
                dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
            end
        end
    end
end

