function Plot
RSDAH = importdata('RSDAveAH.mat');
RCGAH = importdata('RCGAveAH.mat');
LRBFGSAH = importdata('LRBFGSAveAH.mat');
BINIAH = importdata('BINIAveAH.mat');
AH = importdata('initial_time.mat');

%% Plot time VS averaged distance
figure(1)
semilogy(1,100)
hold on
set(gca,'FontSize',13)
RSDAHx = RSDAH.RSD_dis_t+AH;
RCGAHx = RCGAH.RCG_dis_t+AH;
LRBFGSAHx = LRBFGSAH.LRBFGS_dis_t+AH;
fRSDAH = stairs(RSDAHx,RSDAH.RSD_dis,'r','linewidth',2);
fRCGAH = stairs(RCGAHx,RCGAH.RCG_dis,'b','linewidth',2);
fLRBFGSAH = stairs(LRBFGSAHx,LRBFGSAH.LRBFGS_dis,'k','linewidth',2);
fBINIAH = stairs(BINIAH.BINI_dis_t+AH, BINIAH.BINI_dis,'g','linewidth',2);
l = legend([fRSDAH,fRCGAH,fLRBFGSAH,fBINIAH],'RSD','RCG','LRBFGS','RL Iteration');
set(l,'Interpreter','Latex');
xlabel('time (s)','FontSize',16,'Interpreter','Latex');
ylabel('$dist(\mu, X_t)$','FontSize',16,'Interpreter','Latex');
r = max(RSDAH.RSD_dis_t(end),RCGAH.RCG_dis_t(end));
r = max(r,BINIAH.BINI_dis_t(end)) + AH;
axis([0 r 1e-16 1e2])



%% Plot iterations VS distance
figure(2)
semilogy(1,100)
hold on
set(gca,'FontSize',13)
fRSDAH = plot(RSDAH.RSD_dis_iter,'r-*','markersize',12,'linewidth',2);
fRCGAH = plot(RCGAH.RCG_dis_iter,'b-square','markersize',12,'linewidth',2);
fLRBFGSAH = plot(LRBFGSAH.LRBFGS_dis_iter,'k-d','markersize',12,'linewidth',2);
fBINIAH = plot(BINIAH.BINI_dis_iter,'g-v','markersize',12,'linewidth',2);
l = legend([fRSDAH,fRCGAH,fLRBFGSAH,fBINIAH],'RSD','RCG','LRBFGS','RL Iteration');
set(l,'Interpreter','Latex');
xlabel('iterations','FontSize',16,'interpreter','latex');
ylabel('$dist(\mu, X_k)$','FontSize',16,'interpreter','latex');

delete('RSDAveAH.mat','RCGAveAH.mat','LRBFGSAveAH.mat','BINIAveAH.mat','initial_time.mat');

end