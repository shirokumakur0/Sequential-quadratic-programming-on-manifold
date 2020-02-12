function fs = SPD(retraction,vector_transport)

% This function return some necessary function handles of the SPD matrix manifold
% 
% INPUT:
% Grassmann manifold : Gr(p, n)
%     params.retraction : 
%         '1', exponential mapping
%         '2', retraction (second order approximation to the exponential mapping)
%     params.vector_transport : the isometric vector transport
%         '1', vector transport parallelization
%         '2', construct vector transport from retraction 2 that satisfies the locking
%         condition 
% 
% OUTPUT:
%     fns.STV(a, v) : return a tangent vector which is a times v.
%     fns.TVaddTV(v1, v2) : return a tangent vector which is equal to v1
%     fns.inpro(x, v1, v2) : return the inner product g_x(v1, v2) of two tangent vectors v1 and v2 on T_x M.
%     fns.Hv(x, H, v) : return tangent vector which is given by that operator H affect on tangent vector v on T_x M.
%     fns.proj(x, eta) : return P_x(eta) by projecting v to tangent space of x.
%     fns.R(x, eta) : return R_x(eta), where R is a retraction, x is an element on the manifold and eta is a tangent vector of T_x M.
%     fns.beta(x, eta) : return \|eta\| / \|TR_eta (eta)\|, where TR is the differentiated retraction.
%     fns.Tranv(x1, d, x2, v) : return a tangent vector on T_x2 M, which is given by vector transport that transport v \in T_x M to T_{R_x(d)} M. Here x2 = R_x(d).
%     fns.invTranv(x1, d, x2, v) : return a tangent vector on T_x1 M, which is given by inverse vector transport that transport v \in T_x2 M to T_x1 M. Here x2 = R_x1(d).
%     fns.TranH(x1, d, x2, H1) : return a operator, H2, that affects on T_x2 M. H2 satisfies that for any v \in T_x2(M), H2(v) = Tran (H1 (Tran^{-1}(v)))
%                             where Tran is the vector transport fns.Tranv(x1, d, x2, v).
%     fns.rank1operator(x, v1, v2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = g_x(v2, v) * v1.
%     fns.operadd(x, H1, H2) : return a operator, H, on T_x M that satisfies that for any v \in T_x M, H(v) = H1(v) + H2(v).
% 



cprintf([0.1,0.7,0.8], sprintf('Symmetric positive-definite matrix manifold\n'))
fs.proj = @proj; 
fs.operadd = @operadd; 
fs.Hv = @Hv; 
fs.STV = @STV; 
fs.TVaddTV = @TVaddTV;
fs.rank1operator = @rank1operator; 
fs.beta = @beta_parallel; 
fs.inpro = @inpro_intrinsic;


if(retraction == 1)
    fs.R = @Exp;
elseif(retraction == 2)
    fs.R = @Retraction;
end


if(vector_transport == 1)
    fs.Tranv = @parallel;
    fs.invTranv = @parallel;
    fs.TranH = @TranH_pal;
elseif(vector_transport == 2)
    if(retraction == 2)
        fs.Tranv = @Tranv_locking;
        fs.invTranv = @invTranv_locking;
        fs.TranH = @TranH_locking;
        fs.beta = @beta_locking;
    else
        fprintf('convergence is not guarantee');
        return;
    end
end
end


% intrinsic representation for inner product
function output = inpro_intrinsic(x, v1, v2)
output = v1.TV' * v2.TV;
end


% exponential mapping
function output = Exp(x,eta) 
feta = full_pd(x,eta);
temp = x.invL * feta.TV * x.invL';
output.U = x.L * expm(temp) * x.L';
output.U = 0.5 * (output.U + output.U');
output.L = chol(output.U,'lower');
output.invL = pinv(output.L);
end


% retraction
function output = Retraction(x, eta) 
    feta = full_pd(x, eta);
    %xinvU = x.invL' * x.invL;
    output.U = x.U + feta.TV + 0.5 * feta.TV * x.invU * feta.TV;
    output.U = 0.5 * (output.U + output.U');
    output.invU = pinv(output.U);
    output.L = chol(output.U,'lower');
end


% vector transport by parallelization
function [v, d, x, y] = parallel(x, d, y, v)%vector transport
end


% isometric vector transport with locking condition
function [output, d, x, y] = Tranv_locking(x, d, y, v) 
    if(~isfield(d, 'nu1') || ~isfield(d, 'nu2'))
        fd = full_pd(x, d);
        eps1 = d;
        Tee = Tranv_R(x, fd, y, fd);
        Tee_intr = vech_pd(y, Tee);
        d.beta = sqrt(inpro_intrinsic(x, d, d) / inpro_intrinsic(y, Tee_intr, Tee_intr));
        eps2.TV = d.beta * Tee_intr.TV;
        d.nu1.TV = 2 * eps1.TV;
        d.nu2.TV = - eps1.TV - eps2.TV;
    end
    nu1 = d.nu1;
    nu2 = d.nu2;
    output.TV = v.TV - 2 * inpro_intrinsic(y, nu1, v) / inpro_intrinsic(y, nu1, nu1) * nu1.TV;
    output.TV = output.TV - 2 * inpro_intrinsic(y, nu2, output) / inpro_intrinsic(y, nu2, nu2) * nu2.TV;
    output.TV = reshape(output.TV, size(d.TV));
end



% inverse of vector transport with locking condition
function [output, d, x, y] = invTranv_locking(x, d, y, v) 
    if(~isfield(d, 'nu1') || ~isfield(d, 'nu2'))
        fd = full_pd(x, d);
        eps1 = d;
        Tee = Tranv_R(x, fd, y, fd);
        Tee_intr = vech_pd(y, Tee);
        d.beta = sqrt(inpro_intrinsic(x, d, d) / inpro_intrinsic(y, Tee_intr, Tee_intr));
        eps2.TV = d.beta * Tee_intr.TV;
        d.nu1.TV = 2 * eps1.TV;
        d.nu2.TV = - eps1.TV - eps2.TV;
    end
    nu1 = d.nu1;
    nu2 = d.nu2;
    v.TV = reshape(v.TV, [], 1);
    output.TV = v.TV - 2 * inpro_intrinsic(y, nu2, v) / inpro_intrinsic(y, nu2, nu2) * nu2.TV;
    output.TV = output.TV - 2 * inpro_intrinsic(y, nu1, output) / inpro_intrinsic(y, nu1, nu1) * nu1.TV;
end


function output = beta_locking(x1, eta, x2)
    output = eta.beta;
end


function output = TranH_locking(x, d, y, H)
    vv1.TV = reshape(d.nu1.TV, [], 1);
    vv2.TV = reshape(d.nu2.TV, [], 1);
    output = H - 2 * (H * vv1.TV) * (vv1.TV' / (vv1.TV' * vv1.TV));
    output = output - 2 * (output * vv2.TV) * (vv2.TV' / (vv2.TV' * vv2.TV));
    output = output - 2 * vv1.TV * ((vv1.TV' * output) / (vv1.TV' * vv1.TV));
    output = output - 2 * vv2.TV * ((vv2.TV' * output) / (vv2.TV' * vv2.TV));
end


% the approximation 
function output = TranH_pal(x, d, y, H)
output = H;
end



function [output, d, x, y] = Tranv_R(x, d, y, v) % Differentiated vector transport
    temp = 0.5 * v.TV * x.invU * d.TV;
    output.TV = v.TV + temp + temp';
end


% rank-one operator
function output = rank1operator(x, v1, v2)
output = v1.TV * v2.TV';
end


function output = beta_parallel(x1, eta, x2)
    output = 1.0;
end



%% Functions

function output = proj(x, eta)
    output = eta;
end


function output = operadd(x, H1, H2)
    output = H1 + H2;
end


function output = Hv(x, H, v)
    output.TV = H * v.TV;
end


function output = STV(x, a, eta)
    output.TV = a * eta.TV;
end


function output = TVaddTV(x, eta1, eta2)
    output.TV = eta1.TV + eta2.TV; 
end


% extrinsic representation
function output = full_pd(x, v) 
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


% intrinsic representation
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
