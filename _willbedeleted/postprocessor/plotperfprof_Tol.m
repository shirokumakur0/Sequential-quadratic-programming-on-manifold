function plotperfprof_Tol
    

    
    ftol = 1.02; % the tolerance ratio of function value
    constrtol = 1; % Max violation of constraint.


    
    % AO
    %filenames = ["with_SQP_zz_AO_Dim10.dat",...
    %    "with_SQP_zz_AO_Dim15.dat", "with_SQP_zz_AO_Dim20.dat","with_SQP_zz_AO_Dim25.dat",...
    %    "with_SQP_zz_AO_Dim30.dat"];
    %titlenames = ["Dimension 10", "Dimension 15",...
    %    "Dimension 20","Dimension 25","Dimension 30"];
    %locs = {'southeast','southeast','southeast','southeast','southeast','southeast'};
    %set_num = 5;
    
    % BC
    %filenames = ["with_SQP_zz_BC_Dim10.dat",...
    %    "with_SQP_zz_BC_Dim15.dat", "with_SQP_zz_BC_Dim20.dat","with_SQP_zz_BC_Dim25.dat",...
    %    "with_SQP_zz_BC_Dim30.dat"];
    %titlenames = ["Dimension 10", "Dimension 15",...
    %    "Dimension 20","Dimension 25","Dimension 30"];
    %locs = {'southeast','southeast','southeast','southeast','southeast','southeast'};
    %set_num = 5;
    
    %NNPCA
    filenames = ["with_SQP_zz_NNPCA_Dim30.dat", "with_SQP_zz_NNPCA_Dim50.dat",...
        "with_SQP_zz_NNPCA_Dim100.dat","with_SQP_zz_NNPCA_Dim150.dat"], "with_SQP_zz_NNPCA_Dim200.dat"] %, ...
        %"with_SQP_zz_NNPCA_Dim25.dat","with_SQP_zz_NNPCA_Dim30.dat","with_SQP_zz_NNPCA_Dim50.dat"];
    titlenames = ["Dimension 30", "Dimension 50", "Dimension 100","Dimension 150",...
        "Dimension 200"]% , "Dimension 25", "Dimension 30", "Dimension 50"];
    locs = {'southeast','southeast','southeast','southeast','southeast','southeast','southeast'};
    set_num = 5;
    
    
    startingsolver = [1 1 1 1 1 1];
    fig = figure;
    
    for plotnum = 1:set_num % 7 NNPCA % 6 BC, AO
    
        filename = filenames{plotnum};

        data = csvread(filename);
        [nrow, ~] = size(data);

        T = zeros(nrow, 6);
        for ii = 1 : nrow
            extable = data(ii, 1 : 18);
            extable = extable';
            extable = reshape(extable, [3, 6]);
            extable(extable == 0) = eps;
            [T(ii, :), ~,~,~] = timeplotprof(extable, ftol, constrtol);
        end
        %subplot(4,4,plotnum);  % RCs
        %subplot(3,2,plotnum);  % BC, AO
        subplot(3,2,plotnum); %NNPCA
        perf(T(:,startingsolver(plotnum):6), 1, plotnum);
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
        r(find(isnan(r))) = 2*(max_ratio);
        r = sort(r);
        
        disp(r)
        
        %if ns == 7
        %    r = circshift(r, [0,-1]);
        %end
        %if ns == 7
        %    for s = 1: ns
        %        [xs, ys] = stairs(r(:,s), [1 : np]/np);
        %        plot(xs, ys, 'LineStyle', lines{s});
        %        hold on;
        %    end
        %else
            for s = 1: ns
                [xs, ys] = stairs(r(:,s), [1 : np]/np);
                plot(xs, ys, 'LineStyle', lines{s});
                hold on;
            end            
        %end
        
            
        axis([ -0.1 1*(max_ratio)+0.01 0 1 ]);
        %axis([ 0 2*(max_ratio)+0.01 0 1 ]);
        
        
        %if ns == 7
        %    legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)', 'fmincon_interior_point', 'fmincon_SQP','Riemannian SQP','REPMSD'},'Location',locs{plotnum},'Interpreter','latex');
        %else
            legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)', 'fmincon_interior_point','fmincon_SQP','Riemannian SQP'},'Location',locs{plotnum},'Interpreter','latex');
        %end
        ylabel('Performance Profile');
        xlabel('$$\log_2\tau$','Interpreter','latex');
        title(titlenames(plotnum));
    end

    %response = fig2plotly(fig, 'filename', 'matlab-basic-line');
    %plotly_url = response.url;

end
