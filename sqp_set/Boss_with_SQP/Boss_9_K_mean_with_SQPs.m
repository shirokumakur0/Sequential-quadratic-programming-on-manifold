%%
%---------------------------------K means------------------------------------
close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

data_table = {'iris.csv',  4, 3;
              'wine.csv', 13, 2;
              'sonar.csv', 60, 2;
              'ecoli.csv', 5, 2;
              'pima.csv',  8, 2;
              'vehicle.csv', 18, 4;
              'cloud.csv', 10, 4};


n_repeat = 4;   % Number of repeat experiment

for repeat = 1: n_repeat
    
    for dataset = 1 : length(data_table)
        filename = data_table{dataset, 1};
        rank = data_table{dataset, 3};
        
        %_______Set up data______
        data = csvread(filename);
        D = data*data.';

        %________Experiment_____
        options.maxOuterIter = 5000;
        options.maxtime = 3600;
        options.minstepsize = 1e-5;
        options.verbosity = 1;
        %________Setting________
        setting.filename = filename;
        setting.rank = rank;
        setting.maxOuterIter = options.maxOuterIter;
        setting.maxtime = options.maxtime;
        setting.minstepsize = options.minstepsize;
        setting.verbosity = options.verbosity;
        
        % Only do mini-sum-max for low dimensional data
        if dataset <= 1
            specifier.ind = ones(7,1);
        else
            specifier.ind = [0, 1, 1, 1, 1, 1, 1];
        end
        
        result = clientconstraint_stiefel_Kmeans_with_SQP(D, rank, options, specifier, setting);
        result = result(:);
        param = [dataset; repeat];
        outputdata = [result; param]';
        
        filename = sprintf('with_SQP_zz_KM.dat');
        dlmwrite(filename, outputdata, 'delimiter', ',', 'precision', 16, '-append');
    end
    
end
    
    