% Check if exceeding the rank constraints
% For RC_nn
row_dim = [2,3,4,5,6];
for rdim = row_dim
    
    if rdim == 6
        tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance
    else
        tolKKTrespowerset = [2, 3, 4, 5, 6, 7, 8]; % 1e-* tolerance
    end

    if rdim == 2 || rdim == 3
        col_ratio = [1.5 ,2];
    else
        col_ratio = [1.5, 1.75 ,2];
    end
    
    for tolKKTres = tolKKTrespowerset

        constrtol = 10^(-tolKKTres) + 1e-10; % Max violation of constraint.

        for cratio = col_ratio
            cdim = ceil(cratio * rdim);
            %filename = sprintf("with_SQP_zz_RC_lc_RDim%dCDim%dTol%d.dat", rdim, cdim, tolKKTres);
            filename = sprintf("with_SQP_zz_RC_nn_RDim%dCDim%dTol%d.dat", rdim, cdim, tolKKTres);

            data = csvread(filename);
            [nrow, ~] = size(data);
            
            for ii = 1:nrow
                    extable = data(ii, 1 : 12);
                    extable = extable';
                    extable = reshape(extable, [3, 4]);
                for solver = 1:4
                    if isnan(extable(1,solver))
                        fprintf("rank constraints is broken: Rdim%dCdim%dTol%dSolver%d\n",rdim,cdim,tolKKTres,solver);
                    end
                end
            end
        end
    end
end

% Check if exceeding the rank constraints
% For RC_lc
row_dim = [2,3,4];
tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance

for rdim = row_dim

    if rdim == 2 || rdim == 3
        col_ratio = [1.5 ,2];
    else
        col_ratio = [1.5, 1.75 ,2];
    end
    
    for tolKKTres = tolKKTrespowerset

        for cratio = col_ratio
            cdim = ceil(cratio * rdim);
            filename = sprintf("with_SQP_zz_RC_lc_RDim%dCDim%dTol%d.dat", rdim, cdim, tolKKTres);

            data = csvread(filename);
            [nrow, ~] = size(data);
            
            for ii = 1:nrow
                    extable = data(ii, 1 : 12);
                    extable = extable';
                    extable = reshape(extable, [3, 4]);
                for solver = 1:4
                    if isnan(extable(1,solver))
                        fprintf("rank constraints is broken: Rdim%dCdim%dTol%dSolver%d\n",rdim,cdim,tolKKTres,solver);
                    end
                end
            end
        end
    end
end
