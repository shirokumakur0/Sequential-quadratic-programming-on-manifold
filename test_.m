clear;
m = 5;
n = 10;
k = 2;
problem.M = fixedrankembeddedfactory(m, n, k);

xCur = problem.M.rand();

mat = xCur.U * xCur.S * xCur.V.';

