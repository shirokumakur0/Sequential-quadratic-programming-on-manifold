% --- RFW for fast Karcher mean --- %
% author: M. Weber
% input: PSD matrices A{1}, A{2},..., A{n}; x0 (initialization); mit (max number of iterations)
% returns: x (Karcher mean), loss (loss at each iteration), time (cpu time at each iteration), grad (gradient norm at each iteration)

function [x, loss, time, grad]=FWR(A,x0,mit)

   toc_rec=0;
    
   tic;

   n = numel(A);
   d = size(A{1},1);

   [hm Ai] = harmonicMean(A);
   x  = x0;
   am = arithmeticMean(A);
   
   toc_rec=toc_rec+toc;
   
   for k=1:mit
      tic;
      x=0.5*(x+x'); % force symmetric
      s = gradient(x,Ai);
      zk = FWR_dir(x,s,hm,am);
      al = 0.5/k; %2/(k+2);
      xh=x^-0.5;
      xh2=x^0.5;
      x  = xh2*(xh*zk*xh)^al *xh2;
      
      toc_rec=toc_rec+toc;
      
      % record results
      loss(k)=gmobj(x,A);
      time(k)=toc_rec;
      grad(k)=norm(s,'fro');
      
   end
end
