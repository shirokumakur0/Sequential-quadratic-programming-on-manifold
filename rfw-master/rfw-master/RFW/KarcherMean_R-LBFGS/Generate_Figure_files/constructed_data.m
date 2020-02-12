%% data construction

function  [Kmean, As, cond_Kmean, cond_As] = constructed_data(k,n,f)

% generate true Karcher mean
Kmean = eye(n);
Kmeanh = Kmean^(0.5);
Kmeaninvh = Kmean^(-0.5);
cond_Kmean = cond(Kmean);

% generate directions such that their sum is 0
eta = cell(k,1);
W = cell(k,1);


for i = 1 : k
    [O,~] = qr(randn(n));
    m = floor(n/3);
    D{i} = [rand(1, m)+1,(rand(1, n-m)+1)*10^(-f)];
    D{i} = diag(D{i});
    D{i} = D{i}/norm(D{i},2);
    W{i} = O * D{i} * O';
    eta{i} = Kmeanh * logm(Kmeaninvh*W{i}*Kmeaninvh) * Kmeanh;
    eta{i} = 0.5 * (eta{i} + eta{i}');
end
eta = force_eta(eta, k);


% construct the data point
As = cell(k,1);
for i = 1 : k
    temp = Kmeaninvh*eta{i}*Kmeaninvh;
    [V,D] = eig(temp);
    temp1 = V * diag(exp(diag(D))) * V';
    As{i} = Kmeanh*temp1*Kmeanh;
    As{i} = (As{i} + As{i}')/2;
    cond_As(i) = cond(As{i});
end

disp('conditioning of A')
for i = 1 : k
    disp(cond_As(i))
end

end



function eta = force_eta(eta, k)

if(k == 3)
        eta{3} = -eta{1} - eta{2};
end

if(k == 30) 
        for i = 1 : 3
            j = 10*(i-1);
            eta{j+1} = -eta{j+2} - 0.5*eta{j+3};
            eta{j+5} = -eta{j+4} - 0.5*eta{j+3};
            eta{j+6} = -eta{j+7} - 0.5*eta{j+8};
            eta{j+10} = -eta{j+9} - 0.5*eta{j+8};
        end    
end

if(k == 100) 
        for i = 1 : 10
            j = 10*(i-1);
            eta{j+1} = -eta{j+2} - 0.5*eta{j+3};
            eta{j+5} = -eta{j+4} - 0.5*eta{j+3};
            eta{j+6} = -eta{j+7} - 0.5*eta{j+8};
            eta{j+10} = -eta{j+9} - 0.5*eta{j+8};
        end    
end
end


