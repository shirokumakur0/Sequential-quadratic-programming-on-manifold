% --- Computes geometric mean of two matrices --- %
% author: S. Sra
% input: PSD matrices a, b
% returns: geometric mean

function m = gm(a, b)

  ah = sqrtpsd(a);
  iah = inv(ah);
  t  = sqrtpsd(iah*b*iah);
  m  = ah*t*ah;
end
