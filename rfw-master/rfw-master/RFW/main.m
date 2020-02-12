% -------------------------------- %
%   RFW for fast geometric mean    %
% -------------------------------- %
% author: M. Weber

clear all;
close all;
rng(5);

cd('./manopt');
importmanopt;
cd('..');

%% Parameter
n=50; % dim
m=100; % #size
K=20; % #max iterations

%% Generate collection of PSD
A=genPosdef(n,m);
cond_no = cond_numbers(A);

% generate weights
w=rand(1,m);
s=sum(w);
w=w./s;

% means
am=arithmeticMean(A);
[hm Ai] =harmonicMean(A);

% initialization
means{1}=am;
means{2}=hm;
x0 = hm;

% record running time
toc_zhang=0;
toc_E=0;
toc_R=0;
toc_mm=0;

%% Manopt
%disp('************* Running MANOPT CG *************');
%man1=manoptGM(A,x0,K,1); % use CG
%disp('************* Running MANOPT SD *************');
%man2=manoptGM(A,x0,K,2); % use SD
disp('************* Running BB method *************');
man4=manoptGM(A,x0,K,5); % use BB

%% LBFGS
disp('************* Running LBFGS *************');
ImportGROPT
[X_lb, F_lb, G_lb, T_lb, timecost_lb, iter_lb]= Karcher_mean(A,3,0.0001,x0);

for k=1:length(T_lb)
    out_lb(k)=gmobj(X_lb{k}.U,A);
end

%% MM
[mm,it,param,out_mm,time_mm]=karcher(A,K,x0);

%% Zhang
[meannew iter time_zhang allmeannew_zhang]=find_mean_mm(A,x0,K);

for k=1:length(time_zhang)
    out_zhang(k)=gmobj(allmeannew_zhang{k},A);
end

%% EFW
[x, it, tm, output_E, time_E, grad_E]=FWE(A,x0,K);

%% RFW
[x, it, tm, output_R, time_R, grad_R]=FWR(A,x0,K);


