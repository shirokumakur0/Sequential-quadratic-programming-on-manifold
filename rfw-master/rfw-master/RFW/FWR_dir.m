% --- Solves RFW log-linear oracle --- %
% author: M. Weber
% input: x (current value), s (gradient), hm (harmonic mean), am (arithmetic mean)
% returns: RFW step direction

function z = FWR_dir(x,s,hm,am)
   
   xh=x^-0.5;
   xh2=x^0.5;
  
   [q,lam] = eig(xh2*s*xh2); 
   d=diag(lam);
   d(d>=0)=0;
   d(d<0)=1;
   
   hm1=q'*xh*hm*xh*q;
   am1=q'*xh*am*xh*q;
   p1=chol(real(am1-hm1));
       
   z1=p1'*diag(d)*p1+hm1;
   z=xh2*q*z1*q'*xh2;
   
end
