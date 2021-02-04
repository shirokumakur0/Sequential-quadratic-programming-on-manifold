% --- Computes smallest eigenvalue in set of PSD matrices --- %
% author: M. Weber
% input: PSD matrices A{1}, ..., A{m}; n (dimension); m (number of matrices)
% returns: smallest eigenvalue

function lam = smallest_eval(A,n,m)
    [v,d]= eigs(A{1},n);
    k=length(d);
    
    lam = d(k,k);
    
    for i=2:m
       [v,d]= eigs(A{i},n);
       k=length(d);
       %d(k,k)
       if(d(k,k)<lam)
           lam = d(k,k);
       end
    end

end
