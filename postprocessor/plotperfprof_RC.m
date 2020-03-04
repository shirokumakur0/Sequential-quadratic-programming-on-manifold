function plotperfprof
    

    
    ftol = 1.02; % the tolerance ratio of function value
    constrtol = 5e-4; % Max violation of constraint.
    %constrtol = 10;
    
    % RC_nn
    %filenames = ["with_SQP_zz_RC_nn_RDim2CDim3.dat",...
    %    "with_SQP_zz_RC_nn_RDim2CDim4.dat", "with_SQP_zz_RC_nn_RDim2CDim5.dat",...
    %    "with_SQP_zz_RC_nn_RDim3CDim5.dat","with_SQP_zz_RC_nn_RDim3CDim6.dat",...
    %    "with_SQP_zz_RC_nn_RDim3CDim8.dat","with_SQP_zz_RC_nn_RDim4CDim6.dat",...
    %    "with_SQP_zz_RC_nn_RDim4CDim8.dat","with_SQP_zz_RC_nn_RDim4CDim10.dat",...
    %    "with_SQP_zz_RC_nn_RDim5CDim8.dat","with_SQP_zz_RC_nn_RDim5CDim10.dat",...
    %    "with_SQP_zz_RC_nn_RDim5CDim13.dat","with_SQP_zz_RC_nn_RDim7CDim11.dat",...
    %    "with_SQP_zz_RC_nn_RDim7CDim14.dat","with_SQP_zz_RC_nn_RDim7CDim18.dat"];
    %titlenames = ["Dimension Row:2 Col:3", "Dimension Row:2 Col:4", "Dimension Row:2 Col:5","Dimension Row:3 Col:5",...
    %    "Dimension Row:3 Col:6", "Dimension Row:3 Col:8", "Dimension Row:4 Col:6",...
    %    "Dimension Row:4 Col:8", "Dimension Row:4 Col:10", "Dimension Row:5 Col:8",...
    %    "Dimension Row:5 Col:10", "Dimension Row:5 Col:13", "Dimension Row:7 Col:11",...
    %    "Dimension Row:7 Col:14", "Dimension Row:7 Col:18"];
    %locs = {'southeast','southeast','southeast','southeast','southeast',...
    %    'southeast','southeast','southeast','southeast','southeast',...
    %    'southeast','southeast','southeast','southeast','southeast'};
    %set_num = 15;

    % RC_nnlc
    %filenames = ["with_SQP_zz_RC_nnlc_RDim2CDim3.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim2CDim4.dat", "with_SQP_zz_RC_nnlc_RDim2CDim5.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim3CDim5.dat","with_SQP_zz_RC_nnlc_RDim3CDim6.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim3CDim8.dat","with_SQP_zz_RC_nnlc_RDim4CDim6.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim4CDim8.dat","with_SQP_zz_RC_nnlc_RDim4CDim10.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim5CDim8.dat","with_SQP_zz_RC_nnlc_RDim5CDim10.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim5CDim13.dat","with_SQP_zz_RC_nnlc_RDim7CDim11.dat",...
    %    "with_SQP_zz_RC_nnlc_RDim7CDim14.dat","with_SQP_zz_RC_nnlc_RDim7CDim18.dat"];
    %titlenames = ["Dimension Row:2 Col:3", "Dimension Row:2 Col:4", "Dimension Row:2 Col:5","Dimension Row:3 Col:5",...
    %    "Dimension Row:3 Col:6", "Dimension Row:3 Col:8", "Dimension Row:4 Col:6",...
    %    "Dimension Row:4 Col:8", "Dimension Row:4 Col:10", "Dimension Row:5 Col:8",...
    %    "Dimension Row:5 Col:10", "Dimension Row:5 Col:13", "Dimension Row:7 Col:11",...
    %    "Dimension Row:7 Col:14", "Dimension Row:7 Col:18"];
    %locs = {'southeast','southeast','southeast','southeast','southeast',...
    %    'southeast','southeast','southeast','southeast','southeast',...
    %    'southeast','southeast','southeast','southeast','southeast'};
    %set_num = 15;
    
    % RC_nnlc
    filenames = ["with_SQP_zz_RC_lc_RDim2CDim3.dat",...
        "with_SQP_zz_RC_lc_RDim2CDim4.dat", "with_SQP_zz_RC_lc_RDim2CDim5.dat",...
        "with_SQP_zz_RC_lc_RDim3CDim5.dat","with_SQP_zz_RC_lc_RDim3CDim6.dat",...
        "with_SQP_zz_RC_lc_RDim3CDim8.dat","with_SQP_zz_RC_lc_RDim4CDim6.dat",...
        "with_SQP_zz_RC_lc_RDim4CDim8.dat","with_SQP_zz_RC_lc_RDim4CDim10.dat",...
        "with_SQP_zz_RC_lc_RDim5CDim8.dat","with_SQP_zz_RC_lc_RDim5CDim10.dat",...
        "with_SQP_zz_RC_lc_RDim5CDim13.dat","with_SQP_zz_RC_lc_RDim7CDim11.dat",...
        "with_SQP_zz_RC_lc_RDim7CDim14.dat","with_SQP_zz_RC_lc_RDim7CDim18.dat"];
    titlenames = ["Dimension Row:2 Col:3", "Dimension Row:2 Col:4", "Dimension Row:2 Col:5","Dimension Row:3 Col:5",...
        "Dimension Row:3 Col:6", "Dimension Row:3 Col:8", "Dimension Row:4 Col:6",...
        "Dimension Row:4 Col:8", "Dimension Row:4 Col:10", "Dimension Row:5 Col:8",...
        "Dimension Row:5 Col:10", "Dimension Row:5 Col:13", "Dimension Row:7 Col:11",...
        "Dimension Row:7 Col:14", "Dimension Row:7 Col:18"];
    locs = {'southeast','southeast','southeast','southeast','southeast',...
        'southeast','southeast','southeast','southeast','southeast',...
        'southeast','southeast','southeast','southeast','southeast'};
    set_num = 15;

    
    
    startingsolver = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    fig = figure;
    
    for plotnum = 1:set_num % 7 NNPCA % 6 BC, AO
    
        filename = filenames{plotnum};

        data = csvread(filename);
        [nrow, ~] = size(data);

        T = zeros(nrow, 4);
        for ii = 1 : nrow
            extable = data(ii, 1 : 21);
            extable = extable';
            extable = reshape(extable, [3, 7]);
            extable(extable == 0) = eps;
            extable = [extable(:,2:4), extable(:,7)] % added
            [T(ii, :), ~,~,~] = timeplotprof(extable, ftol, constrtol);
        end
                
        subplot(4,4,plotnum);  % RCs
        %subplot(3,2,plotnum);  % BC, AO
        %subplot(4,2,plotnum); %NNPCA
        perf(T(:,startingsolver(plotnum):4), 1, plotnum);
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
            for s = 1: ns
                [xs, ys] = stairs(r(:,s), [1 : np]/np);
                plot(xs, ys, 'LineStyle', lines{s});
                hold on;
            end
        %else
        %    for s = 1: ns
        %        [xs, ys] = stairs(r(:,s), [1 : np]/np);
        %        plot(xs, ys, 'LineStyle', lines{s});
        %        hold on;
        %    end            
        %end
        
            
        axis([ -0.1 1*(max_ratio)+0.01 0 1 ]);
        %axis([ 0 2*(max_ratio)+0.01 0 1 ]);
        
        
        %if ns == 7
        %    legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)', 'fmincon_interior_point', 'fmincon_SQP','Riemannian SQP','REPMSD'},'Location',locs{plotnum},'Interpreter','latex');
        %else
            legend({'RALM','REPMS($Q^{lqh}$)', 'REPMS($Q^{lse}$)','Riemannian SQP'},'Location',locs{plotnum},'Interpreter','latex');
        %end
        ylabel('Performance Profile');
        xlabel('$$\log_2\tau$','Interpreter','latex');
        title(titlenames(plotnum));
    end

    %response = fig2plotly(fig, 'filename', 'matlab-basic-line');
    %plotly_url = response.url;

end
