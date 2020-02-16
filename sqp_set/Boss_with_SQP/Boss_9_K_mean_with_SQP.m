%%
%---------------------------------K means------------------------------------
close all; clear all; clc;
specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

% data_table = {'iris.csv',  4, 3;
%               'wine.csv', 13, 2;
%               'sonar.csv', 60, 2;
%               'ecoli.csv', 5, 2;
%               'pima.csv',  8, 2;
%               'vehicle.csv', 18, 4;
%               'cloud.csv', 10, 4};
data_table = { % data number, data dimesion, clusters
               10,   5,  3;
               25,   8,  3;
               25,   6,  4;
               30,   8,  4;
               35,   7,  5;
               40,   7,  5;
               50,  10,  4;
               60,  10,  5;
               80,  10,  5;
               100,  50,  4;
               200,  75,  5;
               500, 100,  5};
          
[row_data_table,~] = size(data_table);
n_repeat = 1;   % Number of repeat experiment

for repeat = 1: n_repeat
    
    for dataset = 1 : row_data_table
        data_num = data_table{dataset, 1};
        data_dim = data_table{dataset, 2};
        rank = data_table{dataset, 3};
        
        %_______Set up data______
        data = 5 * rand(data_num,data_dim);
        D = data*data.';
        D = 0.5 * (D + D.');
        
        %________Experiment_____
        options.maxOuterIter = 5000;
        options.maxtime = 3600;
        options.minstepsize = 1e-6;
        options.mineigval_correction = 100;
        options.verbosity = 1;
        
        %________Setting________
        setting.data_num = data_num;
        setting.data_dim = data_dim;
        setting.rank = rank;
        setting.repeat = repeat;
        setting.maxOuterIter = options.maxOuterIter;
        setting.maxtime = options.maxtime;
        setting.minstepsize = options.minstepsize;
        setting.verbosity = options.verbosity;
        setting.data = data;
        
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
    
    