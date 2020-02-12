% --- EFW for fast Karcher mean --- %
% author: S. Sra, M. Weber
% input: PSD matrices A{1}, A{2},..., A{n}; x0 (initialization); mit (max number of iterations)
% returns: x (Karcher mean), loss (loss at each iteration), time (cpu time at each iteration), grad (gradient norm at each iteration)

function [x, loss, time, grad]=FWE(A,x0,mit)

   toc_rec=0;
   tic;

   n = numel(A);
   d = size(A{1},1);

   [hm Ai] = harmonicMean(A);
   x  = x0;
   am = arithmeticMean(A);
   p = chol(am-hm);
   
   toc_rec=toc_rec+toc;
   
   for k=1:mit
      tic;
      x=0.5*(x+x'); % force symmetric
      
      s = gradient(x,Ai);
      zk = FWE_dir(s,hm,am,p);
      al = 2/(k+2);
      x  = x + al*(zk-x);
      
      toc_rec=toc_rec+toc;
      
      % record results
      loss(k)=gmobj(x,A);
      time(k)=toc_rec;
      grad(k)=norm(s,'fro');
      
   end
end
