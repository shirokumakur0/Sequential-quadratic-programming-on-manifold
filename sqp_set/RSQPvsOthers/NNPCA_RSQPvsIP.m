%%
%---------------------------------PCA------------------------------------

close all; clc; clear all;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

% Set for performance profile RSQP vs IP
dim_set = [150];
snrset = [0.5];  % Signal Strength
deltaset = [0.7];  % Sparsity 
tolKKTrespowerset = [20]; % 1e-* tolerance

rank = 1;                                  % Rank of BM Relaxation. 1 if we don't.
n_repeat = 1;                              % Number of repeat experiment

for repeat = 1: n_repeat
    
    for dim = dim_set
        
        for snr = snrset
            
            for delta = deltaset
                
                for tolKKTres = tolKKTrespowerset

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
                    options.maxOuterIter = 100000; % for Riemannian methods & fmincon
                    options.maxiter = options.maxOuterIter;  % for RSQP
                    options.maxtime = 600;
                    options.tolKKTres = 10^(-tolKKTres);
                    options.outerverbosity = 1;
                    options.verbosity = options.outerverbosity;
                    %________Setting________
                    setting.repeat = repeat;
                    setting.dim = dim;
                    setting.snr = snr;
                    setting.delta = delta;
                    setting.rank = rank;
                    setting.tolKKTres = tolKKTres;
                    setting.maxOuterIter = options.maxOuterIter;
                    setting.maxtime = options.maxtime;
                    setting.Z = Z;

                    specifier.ind = [1, 1, 1, 1, 1, 1];

                    result = clientconstraint_sphere_nonnegativePCA_with_SQP(X, rank, options, specifier, setting);
                    result = result(:);
                    param = [dim; snr; delta; repeat; tolKKTres];
                    outputdata = [result; param]';

                    % PP according to dimension and tolKKTres
                    filename = sprintf('with_SQP_zz_NNPCA_Dim%dTol%d.dat', dim, tolKKTres);
                    dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                    % PP according to tolKKTres
                    %filename = sprintf('with_SQP_zz_NNPCA_Tol%d.dat', tolKKTres);
                    %dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
                end
            end
        end
    end
end

