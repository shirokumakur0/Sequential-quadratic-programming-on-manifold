% --- Computes Wasserstein distance from x to set of PSD matrices --- %
% author: S. Sra
% input: PSD matrices A{1}, A{2},..., A{n}; PSD matrix x
% returns: Sum of Wasserstein distances from x to all A_i \in A

function o = wmobj(x,A)
   o = 0;
   x = (x+x')/2;
   for i=1:numel(A)
      o = o + bures(x,A{i})^2;
   end
 end
 
 
