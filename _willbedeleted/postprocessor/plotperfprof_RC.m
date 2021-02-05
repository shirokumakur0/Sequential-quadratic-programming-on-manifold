function plotperfprof_RC
    
    row_dim = [2];
    ftol = 1.02; % the tolerance ratio of function value (Won't be used, just a heritage)
    
    fig = figure;
    startingsolver = [1 1 1 1];
    
    plotnum = 1;
    tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance

    for rdim = row_dim
        %if rdim == 7
        %    tolKKTrespowerset = [2]; % 1e-* tolerance
        %elseif (rdim == 6) || (rdim == 5)
        %    tolKKTrespowerset = [2]; % 1e-* tolerance
        %else
        %    tolKKTrespowerset = [2, 3, 4, 5, 6, 7]; % 1e-* tolerance
        %end
      
        if (rdim == 2) || (rdim == 3)
            col_ratio = [1.5, 2, 3];
        else
            col_ratio = [1.5, 1.75, 2, 3];
        end
        
        for tolKKTres = tolKKTrespowerset
            
            constrtol = 10^(-tolKKTres); % Max violation of constraint.

            for cratio = col_ratio
                cdim = ceil(cratio * rdim);
                %filename = sprintf("with_SQP_zz_RC_lc_RDim%dCDim%dTol%d.dat", rdim, cdim, tolKKTres);
                filename = sprintf("with_SQP_zz_RC_nn_RDim%dCDim%dTol%d.dat", rdim, cdim, tolKKTres);
                titlename = sprintf("Dimension Row:%d Col:%d Residual: 1e-%d", rdim, cdim, tolKKTres);
                locs = sprintf('southeast');
            
                data = csvread(filename);
                [nrow, ~] = size(data);

                T = zeros(nrow, 4);
                for ii = 1 : nrow
                    extable = data(ii, 1 : 12);
                    extable = extable';
                    extable = reshape(extable, [3, 4]);
                    extable(extable == 0) = eps;
                    [T(ii, :), ~,~,~] = timeplotprof(extable, ftol, constrtol);
                end

                
                subplot(6,3,plotnum); % nn rdim == 6
                %subplot(6,2,plotnum); % nn rdim == 3|4|5
                %subplot(7,2,plotnum); % nn rdim == 2
                perf(T, 1, plotnum);
                
                plotnum = plotnum + 1;
            end
        end
    end
    
    function perf(T, logplot, plotnum)
        if (nargin< 2) 
            logplot = 0; 
        end
        
        colors = ['m' 'b' 'r' 'g' 'r' 'k' 'y'];
        co = [0 0 1;
              0 0.5 0;
              1 0 0;
              0 0.75 0.75;
              0.75 0 0.75;
              0.75 0.75 0;
              0.25 0.25 0.25];
        
        lines = {'--' '-' '-.' '-' '--' '--', '--'};
       
        markers = [' ' ' ' ' ' ' ' ' ' '.' '.' ' '];
       
        [np,ns] = size(T);
        minperf = min(T, [], 2);
        r = zeros(np, ns);
        for p = 1: np
            r(p,:) = T(p,:)/minperf(p);
        end
        if (logplot) 
            r = log2(r); 
        end
        
        disp(r)

        max_ratio = max(max(r));
        r(find(isnan(r))) = 1e+50;
        r = sort(r);
        
        disp(r)
        
        for s = 1: ns
            [xs, ys] = stairs(r(:,s), [1 : np]/np);
            plot(xs, ys, 'LineStyle', lines{s});
            hold on;
        end
             
        %axis([ -0.1 1*(max_ratio)+0.01 0 1 ]);
        axis([ -0.001 2*(max_ratio)+0.01 0 1 ]);
       
        legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)','Riemannian SQP'},'Location',locs,'Interpreter','latex');
        ylabel('Performance Profile');
        xlabel('$$\log_2\tau$','Interpreter','latex');
        title(titlename);
    end

    %response = fig2plotly(fig, 'filename', 'matlab-basic-line');
    %plotly_url = response.url;

end
