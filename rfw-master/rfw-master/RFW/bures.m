% --- Computes Bures distance --- %
% author: S. Sra
% input: PSD matrices A, B
% returns: d (Bures distance of A and B), t (cpu time)

function [d t] = bures(A,B)

   t=tic;
   if isvector(A)
      A=diag(A); B=diag(B);
   end
   
   ah=sqrtpsd(A);
   d = trace(A+B)-2*trace(sqrtpsd(ah*B*ah));
   d = sqrt(d);
   t=toc(t);

end
