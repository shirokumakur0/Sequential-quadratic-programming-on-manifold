% --- Computes square root of PSD matrix --- %
% author: S. Sra
% input: PSD matrix X
% returns: square root of X

function L=sqrtpsd(X)
   [u,d]=schur(X);
   d=diag(d);
   d(d<=0)=0;
   L=u*diag(sqrt(d))*u';
end
