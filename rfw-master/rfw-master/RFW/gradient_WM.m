% --- Computes Riemannian gradient for Wasserstein mean --- #
% author: M. Weber
% input: PSD matrices A{1}, A{2},..., A{n}; PSD matrix x
% returns: Riemannian gradient

function s = gradient_WM(x,A)
    id = eye(size(x));
    y = zeros(size(x));
    xi = inv(x);
    n = numel(A);
    for i=1:numel(A)
        y = y + gm(A{i},xi)/n;
    end             
    y = (y+y')/2;
    s = id - y;
end
