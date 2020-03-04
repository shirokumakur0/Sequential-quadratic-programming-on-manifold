% --- Computes Riemannian distance of X to set of PSD matrices A --- %
% author: S. Sra
% input: PSD matrices A{1}, A{2},..., A{n}; PSD matrix x
% returns: Riemannian distance

function o = gmobj(x,A)
   
   o = 0;
   x = (x+x')/2;                
   for i=1:numel(A)
      o = o + riem(x,A{i})^2;
   end
end
