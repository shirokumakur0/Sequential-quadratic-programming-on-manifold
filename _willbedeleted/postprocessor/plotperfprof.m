function plotperfprof    
    
    ftol = 1.02; % the tolerance ratio of function value (Won't be used, just a heritage)
    
    startingsolver = [1 1 1 1 1 1];
    fig = figure;
    
    dim_set = [10, 50, 75, 100];
    tolKKTrespowerset = [2, 4, 6, 8, 10]; % 1e-* tolerance
    
    plotnum = 1;

    for rdim = dim_set      
        for tolKKTres = tolKKTrespowerset
            constrtol = 10^(-tolKKTres); % Max violation of constraint.
            filename = sprintf("with_SQP_zz_BC_Dim%dTol%d.dat", rdim, tolKKTres);
            titlename = sprintf("Dimension: %d Residual: 1e-%d", rdim, tolKKTres);
            locs = sprintf('southeast');

            data = csvread(filename);
            [nrow, ~] = size(data);

            T = zeros(nrow, 5);
            for ii = 1 : nrow
                extable = data(ii, 1 : 18);
                extable(10:12) = [];
                extable = extable';
                extable = reshape(extable, [3, 5]);
                extable(extable == 0) = eps;
                [T(ii, :), ~,~,~] = timeplotprof(extable, ftol, constrtol);
            end

            subplot(4,5,plotnum); % nn rdim == 3|4|5
            perf(T, 1, plotnum);

            plotnum = plotnum + 1;
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
        
        
        legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)','fmincon_SQP','Riemannian SQP'},'Location',locs,'Interpreter','latex');
        ylabel('Performance Profile');
        xlabel('$$\log_2\tau$','Interpreter','latex');
        title(titlename);
    end

    %response = fig2plotly(fig, 'filename', 'matlab-basic-line');
    %plotly_url = response.url;

end
