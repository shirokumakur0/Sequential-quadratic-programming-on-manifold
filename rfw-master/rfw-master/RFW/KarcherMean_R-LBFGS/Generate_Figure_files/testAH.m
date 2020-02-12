function testAH

Flag = 1;

for times = 1 : 5
    if(Flag == 1)
        Data = importdata(['Kmean1' '.mat']);
        As = Data.As;
        Kmean = Data.Kmean;
        k = length(As);
        x1 = Data.AH;
        
        TOL = 1e-2;
        for method = 1:3
            Karcher_mean(As,method,TOL,x1);
        end
        As{k+1} = x1;
        Bini_Karcher(As{1:k+1});
        Flag = 2;
    end
    
    Data = importdata(['Kmean' num2str(times) '.mat']);
    As = Data.As;
    Kmean = Data.Kmean;
    
    k = length(As);
    n = size(Kmean,2);
    x1 = Data.AH;
    
    TOL = 1e-15;
    method = 1; % RSD
    [X, F, G, T, timecost, iter] = Karcher_mean(As,method,TOL,x1);
    [dis, norm2, norminfty] = Dist(Kmean,X);
    Xcom = X{end}.U;
    save(['RSDAH',int2str(times)],'dis','norm2','norminfty','F','G','T','timecost','iter','Xcom');
    
    method = 2; % RCG
    [X, F, G, T, timecost, iter] = Karcher_mean(As,method,TOL,x1);
    [dis, norm2, norminfty] = Dist(Kmean,X);
    Xcom = X{end}.U;
    save(['RCGAH',int2str(times)],'dis','norm2','norminfty', 'F','G','T','timecost','iter','Xcom');
    
    
    method = 3; % LRBFGS
    [X, F, G, T, timecost, iter] = Karcher_mean(As,method,TOL,x1);
    [dis, norm2, norminfty] = Dist(Kmean,X);
    Xcom = X{end}.U;
    save(['LRBFGSAH',int2str(times)],'dis','norm2','norminfty', 'F','G','T','timecost','iter','Xcom');
    
    
    % Bini's Karcher mean
    As{k+1} = x1;
    [X, T, timecost, iter,theta] = Bini_Karcher(As{1:k+1});
    [dis, norm2, norminfty] = Dist(Kmean,X);
    Xcom = X{end}.U;
    save(['BiniAH',int2str(times)],'dis','norm2','norminfty','T','timecost', 'iter','Xcom');
end

end



%% compute distance between two points
function [dis, norm2, norminfty] = Dist(X0,X)
K = length(X);
X0invh = X0^(-0.5);
Xlog = cell(K);
dis = zeros(1,K);
Error = cell(K);
norm2 = zeros(1,K);
norminfty = zeros(1,K);

for i = 1 : K
    [V,D] = eig(X0invh * X{i}.U * X0invh);
    Xlog{i} = V * diag(log(diag(D))) * V';
    dis(i) = norm(Xlog{i}, 'fro'); % Riemannian distance
    Error{i} = X{i}.U - X0;
    norm2(i) = norm(Error{i},2); % matrix 2 norm
    norminfty(i) = norm(Error{i},inf); % matrix infinity norm
end
end