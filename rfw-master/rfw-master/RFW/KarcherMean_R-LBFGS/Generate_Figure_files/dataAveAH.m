function dataAveAH

%% importdata
times = 5;
for i = 1 : times
RSD{i} = importdata(['RSDAH' num2str(i) '.mat']);
RCG{i} = importdata(['RCGAH' num2str(i) '.mat']);
LRBFGS{i} = importdata(['LRBFGSAH' num2str(i) '.mat']);
BINI{i} = importdata(['BiniAH' num2str(i) '.mat']);

delete(['RSDAH' num2str(i) '.mat'],['RCGAH' num2str(i) '.mat'],...
    ['LRBFGSAH' num2str(i) '.mat'], ['BiniAH' num2str(i) '.mat'],...
    ['Kmean' num2str(i) '.mat']);
end


%% compute the average of distance
[RSD_dis,RSD_dis_t] = ave_distance(RSD,times);
[RCG_dis,RCG_dis_t] = ave_distance(RCG,times);
[LRBFGS_dis,LRBFGS_dis_t] = ave_distance(LRBFGS,times);
[BINI_dis,BINI_dis_t] = ave_distance(BINI,times);


%% compute the average of norm2
[RSD_norm2,RSD_norm2_t] = ave_norm2(RSD,times);
[RCG_norm2,RCG_norm2_t] = ave_norm2(RCG,times);
[LRBFGS_norm2,LRBFGS_norm2_t] = ave_norm2(LRBFGS,times);
[BINI_norm2,BINI_norm2_t] = ave_norm2(BINI,times);


%% compute the average of norminfty
[RSD_norminfty,RSD_norminfty_t] = ave_norminfty(RSD,times);
[RCG_norminfty,RCG_norminfty_t] = ave_norminfty(RCG,times);
[LRBFGS_norminfty,LRBFGS_norminfty_t] = ave_norminfty(LRBFGS,times);
[BINI_norminfty,BINI_norminfty_t] = ave_norminfty(BINI,times);


%% compute the average of gradient norm
RSD_grad = ave_grad(RSD,times);
RCG_grad = ave_grad(RCG,times);
LRBFGS_grad = ave_grad(LRBFGS,times);


%% compute the average of distance along with iterations
RSD_dis_iter = ave_dis_iter(RSD,times);
RCG_dis_iter = ave_dis_iter(RCG,times);
LRBFGS_dis_iter = ave_dis_iter(LRBFGS,times);
BINI_dis_iter = ave_dis_iter(BINI,times);

save('RSDAveAH.mat', 'RSD_dis', 'RSD_norm2','RSD_norminfty', 'RSD_grad',...
    'RSD_dis_t','RSD_norm2_t', 'RSD_norminfty_t','RSD_dis_iter');

save('RCGAveAH.mat', 'RCG_dis', 'RCG_norm2','RCG_norminfty', 'RCG_grad',...
    'RCG_dis_t','RCG_norm2_t','RCG_norminfty','RCG_dis_iter');

save('LRBFGSAveAH.mat', 'LRBFGS_dis', 'LRBFGS_norm2','LRBFGS_norminfty', 'LRBFGS_grad',...
    'LRBFGS_dis_t','LRBFGS_norm2_t','LRBFGS_norminfty_t','LRBFGS_dis_iter');

save('BINIAveAH.mat', 'BINI_dis', 'BINI_norm2','BINI_norminfty',...
    'BINI_dis_t','BINI_norm2_t','BINI_norminfty_t','BINI_dis_iter');

end




%% compute the average of norm
function [AveNorm2, t] = ave_norm2(X,times)
Xtemp = [X{:}];
X_max_time = max([Xtemp.T]);

dt = 0.001;
ngrid = floor(X_max_time/dt);
AveNorm2 = zeros(1,ngrid);
t = 0 : dt : (ngrid-1) * dt;

for i = 1: times
    Norm2 = X{i}.norm2(end);%------------
    Grid{i} = Norm2 * ones(1,ngrid);
    for k = 1 : length(X{i}.T)
        ind(k) = floor(X{i}.T(k)/dt);
    end
    Grid{i}(1:ind(1)) = X{i}.norm2(1);
    for j = 1 : length(X{i}.norm2)-1
        Grid{i}(ind(j)+1 : ind(j+1)) = X{i}.norm2(j);
    end
    AveNorm2 = AveNorm2 + Grid{i};
end
AveNorm2 = AveNorm2/times;
end



%% compute the average of norm
function [AveNorminfty, t] = ave_norminfty(X,times)
Xtemp = [X{:}];
X_max_time = max([Xtemp.T]);

dt = 0.001;
ngrid = floor(X_max_time/dt);
AveNorminfty = zeros(1,ngrid);
t = 0 : dt : (ngrid-1) * dt;

for i = 1: times
    Norminfty = X{i}.norminfty(end);%------------
    Grid{i} = Norminfty * ones(1,ngrid);
    for k = 1 : length(X{i}.T)
        ind(k) = floor(X{i}.T(k)/dt);
    end
    Grid{i}(1:ind(1)) = X{i}.norminfty(1);
    for j = 1 : length(X{i}.norminfty)-1
        Grid{i}(ind(j)+1 : ind(j+1)) = X{i}.norminfty(j);
    end
    AveNorminfty = AveNorminfty + Grid{i};
end
AveNorminfty = AveNorminfty/times;
end

%% compute the average of gradient norm
function X_grad = ave_grad(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_min_iter = min([Xtemp.iter])+1;
X_grad = zeros(1,X_max_iter);
for i = 1: times
    Grad = X{i}.G(end);
    X_grad = X_grad + [X{i}.G Grad*ones(1,length(X_grad) - length(X{i}.G))];
end
X_grad = X_grad/times;
end



%% compute the regular average of distance
function X_dis = ave_dis_iter(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_dis = zeros(1,X_max_iter);
for i = 1: times
    dis = X{i}.dis(end);
    X_dis = X_dis + [X{i}.dis dis*ones(1,length(X_dis) - length(X{i}.dis))];
end
X_dis = X_dis/times;
end


%% compute the average of distance
function [AveDis,t] = ave_distance(X,times)
Xtemp = [X{:}];
X_max_time = max([Xtemp.T]);

dt = 0.0001;
ngrid = floor(X_max_time/dt);
AveDis = zeros(1,ngrid);
t = 0 : dt : (ngrid-1) * dt;
for i = 1: times
    Dis = X{i}.dis(end);%------------
    Grid{i} = Dis * ones(1,ngrid);
    for k = 1 : length(X{i}.T)
        ind(k) = floor(X{i}.T(k)/dt);
    end
    Grid{i}(1:ind(1)) = X{i}.dis(1);
    for j = 1 : length(X{i}.dis)-1
        Grid{i}(ind(j)+1 : ind(j+1)) = X{i}.dis(j);
    end
    AveDis = AveDis + Grid{i};
end
AveDis = AveDis/times;
end

