function [X, F, G, T, timecost, iter] = driver_SPD(k, A, method, fns, params)
% driver for geometric mean
%
% INPUT:
% A: the symmetric positive definite matrices in the cost function
% k: the number of the matrices
% n: the size of the matrices
% method: indication of method
% fns: a struct that contains required function handles
% params: a struct that contains parameters
%
% OUTPUT:
% X : all iterates
% F : function values of all iterates
% G : norm of gradient of all iterates
% T : accumulated computational time for iterations
% timecost : total computational time
% nf : number of function evaluations
% ng : number of gradient evaluations
% nR : number of retraction actions
% nH : number of Hessian actions
% nV : number of vector transport actions without parameters
% nVp : number of vector transport action with parameters

fns.f = @(x)f(x, A);
fns.gf = @(x)gf(x, A);
fns.hessianv = @(x,v)hessianv(x,v,A);


nH = 0;
nV = 0;
nVp = 0;

if(method == 1)
    [X, F, G, timecost, nf, ng, T, nV, nR, iter] = RSD(fns, params);
elseif(method == 2)
    [X, F, G, timecost, nf, ng, T, nV, nR, iter] = RCG(fns, params);
elseif(method == 3)
    [X, F, G, timecost, nf, ng, T, nV, nVp, nR, iter] = LRBFGS(fns, params);
elseif(method == 4)
    [X, F, G, timecost, nf, ng, nH, T, nV, nR, iter] = RBFGS(fns, params);
elseif(method == 5)
    [X, F, G, timecost, nf, ng, nH, T, nR, iter] = RTR_Newton(fns, params);
elseif(method == 6)
    [X, F, G, timecost, nf, ng, nH, T, nV, nR, iter] = RNewton(fns, params);
end
end



function [output, x] = f(x, A)  % cost function
k = length(A.invh);
output = 0;
if(sum(find(isinf(x.U))))
    output = inf;
    return;
end

for i = 1 : k
    S = A.invh{i} * x.L;
    [V, D] = schur(S * S');
    x.log{i} = V * diag(log(diag(D))) * V';
    x.c(i) = D(end,end)/D(1,1);
    output = output + norm(x.log{i}, 'fro')^2;
end
output = output / (2 * k);
end


function [output, x] = gf(x, A) % gradient of f(x)
k = length(A.invh);
if(~isfield(x, 'log'))
    for i = 1:k
        S = A.invh{i} * x.L;
        [V, D] = schur(S * S');
        x.log{i} = V * diag(log(diag(D))) * V';
        x.c(i) = D(end,end)/D(1,1);
    end
end
output.TV = zeros(size(x.U));
for i = 1 : k
    output.TV = output.TV + x.U * A.invh{i} * x.log{i} * A.h{i};
end
output.TV = 0.5 * (output.TV + output.TV');
output.TV = output.TV/k;
output = vech_pd(x, output);
end



function output = hessianv(x, v, A)
k = length(A.inv);
output.TV = zeros(size(x.U));
fv = full_pd(x, v);
for i = 1:k
    temp = fv.TV * A.invh{i} * x.log{i} * A.h{i};
    output.TV = output.TV + temp - temp' + ...
        x.U * Dlog(A.inv{i} * x.U, A.inv{i} * fv.TV);
end
output.TV = 0.5 * (output.TV + output.TV');
output.TV = output.TV/k;
output = vech_pd(x, output);
end



function output = full_pd(x, v) % recover v
if(~isfield(x, 'L'))
    x.L = chol(x.U,'lower');
end
n = size(x.U, 1);
omega = tril(ones(n, n), -1);
indx = find(omega);
omega(indx) = v.TV(1 : 0.5 * n * (n - 1)) / sqrt(2);
output.TV = omega + omega' + diag(v.TV(1 + 0.5 * n * (n - 1) : end));
output.TV = x.L * output.TV * x.L';
end



function output = vech_pd(x, eta)
if(~isfield(x, 'invL'))
    x.L = chol(x.U,'lower');
    x.invL = inv(x.L);
end
n = size(eta.TV, 1);
indx = find(tril(ones(n, n), -1));
eta_new = x.invL * eta.TV * x.invL';
output.TV = zeros(0.5 * n * (n + 1), 1);
output.TV(1 : 0.5 * n * (n - 1)) = sqrt(2) * eta_new(indx);
output.TV(1 + 0.5 * n * (n - 1) : end) = diag(eta_new);
end