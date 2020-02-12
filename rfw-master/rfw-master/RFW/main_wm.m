% -------------------------------- %
%  RFW for fast Wasserstein mean   %
% -------------------------------- %
% author: M. Weber

clear all;
close all;
rng(5);

cd('./manopt');
importmanopt;
cd('..');

%% Parameter
n=80; % dim
m=20; % #size
K=60; % #max iterations

%% Generate collection of PSD
A=genPosdef2(n,m);

% means
am=arithmeticMean(A);
[hm Ai] = harmonicMean(A);

% bounds
lmin = smallest_eval(A,n,m);
lb = lmin*eye(size(am));
%lb = n^(2/3)*lmin^(1/3)*eye(size(am));
ub = am;

% initialization
x1 = A{1}; 
x2 = lb;
x3 = am;

%% RFW
[x1, it1, tm1, output_R1, time_R1, grad_R1]=FWR_WM(A,x1,K);
[x2, it2, tm2, output_R2, time_R2, grad_R2]=FWR_WM(A,x2,K);
[x3, it3, tm3, output_R3, time_R3, grad_R3]=FWR_WM(A,x3,K);

%% Plot results
figure(1)
legendStr = {'A_1','\alpha I','AM'};
loglog(time_R1,output_R1,'-','LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])
hold all
loglog(time_R2,output_R2,'-','LineWidth',3,'Color',[0.6350, 0.0780, 0.1840]	)
loglog(time_R3,output_R3,'-','LineWidth',3,'Color',[0.75, 0.75, 0])
legend(legendStr);
xlabel('time (s)');
ylabel('loss function');
set(gca,'linewidth',2);
set(gca,'fontsize', 20);
ax = gca;
ax.FontSize = 22; 
title('Performance');


