function [w,d] = qpsubprob(S, prtlevel)
% Solver for the quadractic convex program
%
%     min          w'S'Sw = ||S w||^2
%     w 
%
%     subject to   e'w = 1          Sw is a convex combination of gradients
%                   w >= 0          in the columns of S
%  Return w and d = -Sw.
%
%  The QP is solved by a call to quadprog which will invoke either the 
%  Matlab Optimization Toolbox or MOSEK software, depending on which is
%  installed and, if both are, their relative order in the Matlab Path.  
%  Generally, MOSEK is far preferable.  To select it, use "addpath" to
%  add it to the front of the path.  To see which is in use, use "which".
%  Written by M. Overton (overton@cs.nyu.edu), last revised November 2004 
%
[m,N] = size(S);  % N gradients of length m
H = S'*S;
e = ones(1,N);    % constraint e'w = 1
O = zeros(1,N);
w0 = ones(N,1)/N;  % starting point for quadprog 
qp_options = optimset('Display','off','TolX', 1e-12, 'TolFun', 1e-12, 'Algorithm', 'active-set');
% qp_options = optimset('Display','off','TolX', 1e-12, 'TolFun', 1e-12, 'LargeScale', 'off');
w = quadprog(H, O, [], [], e, 1, O, [], w0, qp_options);    
% full call: [w, fval, qpflag, qpoutput, lambda] = quadprog(H, O, [], [], e, 1, O, [], w0, qp_options);
if isempty(w)
    error('w is empty: MOSEK license problem?')
end
d = -S*w;
checktol = 1e-7*max(1, norm(H, inf));
if prtlevel > 0 & (min(w) < -checktol | abs(sum(w) - 1) > checktol)
    fprintf('computed w does not satisfy requirements\n')
end
