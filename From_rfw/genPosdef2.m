% --- Generates ill-conditioned PSD matrices ---%
% author: S. Sra
% input: d (dimension), n (number of matrices)
% returns: Set of PSD matrices

function S = genPosdef2(d,n)
   k=d-3;
   delta=0.01;
   for i=1:n
     u = rand(d,k);
     a = delta*eye(d) + u*u';
     S{i} = a;
   end
end
