% --- Computes harmonic mean --- %
% author: S. Sra
% input: PSD matrices A{1}, A{2},..., A{n}
% output: harmonic mean

function [hm, Ai] = harmonicMean(A)
   Ai = A;
   hm = inv(A{1});
   Ai{1}=hm;
   for i=2:numel(A)
      Ai{i} = inv(A{i});
      hm = hm + Ai{i};
   end
   hm = numel(A)*inv(hm);
end

