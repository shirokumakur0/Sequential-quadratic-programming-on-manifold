%%% Solves EFW linear oracle %%%
% author: S. Sra
% input: s (gradient), hm (harmonic mean), am (arithmetic mean), p (factorization)
% returns: step direction

function x = FWE_dir(s,hm,am,p)
    [q,lam] = eig(p*s*p');
    d=diag(lam);
    d(d>=0)=0;
    d(d<0)=1;
    z = q*diag(d)*q';
    x = hm + p'*z*p;
end
