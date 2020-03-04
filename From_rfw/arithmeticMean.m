% --- Computes arithmetic mean --- %
% author: S. Sra
% input: PSD matrices A{1}, A{2},..., A{n}
% returns: arithmetic mean of A

function am = arithmeticMean(A)
   am = A{1};
   for i=2:numel(A)
      am = am + A{i};
   end
   am = am/numel(A);
end
