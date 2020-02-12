% --- Generates well-conditioned PSD matrices --- %
% author: S. Sra
% input: d (dimension), n (number of matrices)
% returns: Set of PSD matrices

function S = genPosdef(d,n)

    for i=1:n
        a = 3*randn(d); a = a*a';
        S{i} = a;
    end

end
