% --- Computes Riemannian distance --- %
% author: S. Sra
% input: PSD matrices A, B
% returns: d (distance), t (cpu time)

function [d t] = riem(A,B)

    t=tic;
    d = eig(A, B);
    d=norm(log(d));

    t=toc(t);

end
