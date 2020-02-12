% --- Computes Riemannian gradient for Karcher mean --- %
% author: S. Sra
% input: PSD matrices Ai{1}, Ai{2},..., Ai{n}; PSD matrix x
% returns: Riemannian gradient

function s = gradient(x,Ai)
   s = zeros(size(x));
   [u,d]=eig(x);
   xh = u*diag(sqrt(diag(d)))*u';
   xhi = u*diag(1./sqrt(diag(d)))*u';
   for i=1:numel(Ai)
      t = xh*Ai{i}*xh;
      t = (t+t')/2;
      [u,d]=eig(t);
      t = u*diag(log(diag(d)))*u';
      s = s + t;
   end
   s = xhi*s*xhi;
   s = (s+s')/2;
end
