
M = fixedrankembeddedfactory(2, 2, 1);
X = struct('U', [1;0], 'V', [1;0], 'S', 1);
V = struct('Up', [0;1], 'Vp', [0;1], 'M', 1);
entry = @(M) M(1, 1);
mat = @(X) X.U*X.S*X.V';
g = @(t) entry(mat(M.retr(X, V, t)));
g(-1)
ezplot(g, [-2, 2]);