% --- Computes condition numbers --- %
% author: M. Weber
% input: PSD matrices A{1}, A{2},..., A{n}
% returns: condition numbers (n-dim array)

function x=cond_numbers(A)
    x=zeros(numel(A));
    for i=1:numel(A)
        x(i)=cond(A{i}); %2-norm cond number
    end
end
    
