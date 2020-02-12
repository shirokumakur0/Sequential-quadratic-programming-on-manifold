% --- Determines if all matrices in set A are positive definite --- %
% author: S. Sra
% input: PSD matrices A{1}, A{2},..., A{m}
% returns: true (if all matrices in A are PSD), false (else)

function r = isposdef(A)
    m=numel(A);
    l=false(m,1);
    for i=1:m
        l(i)=all(eig(A{i}) > eps);
    end
    if all(l==true)
        r=true;
    else
        r=false;
    end
end
