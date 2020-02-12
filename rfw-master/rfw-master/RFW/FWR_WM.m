function [x, it, tm, output, time, grad]=FWR_WM(A,x0,mit)
% FWR for Karcher mean problem
% M. Weber, 08/2017
% input: PSD matrices A{1}, A{2},..., A{n}; x0 (initialization); mit (max number of iterations)

   toc_rec=0;
    
   tic;
   flag = 0;
   if (nargout > 1)
      flag = 1;
      it   = nan*ones(mit,1);
      tm   = nan*ones(mit,1);
   end

   n = numel(A);
   d = size(A{1},1);
   
   % bounds
   am = arithmeticMean(A);
   lmin = smallest_eval(A,d,n);
   lb = lmin*eye(size(am));
   ub = am;
   x  = x0;
   p = chol(real(ub-lb)); 
   
   toc_rec=toc_rec+toc;
   
   for k=1:mit
      tic;
      %x=0.5*(x+x'); % force symmetric
      if (flag) 
         it(k)=wmobj(x,A); 
         tm(k)=toc;
      end
      s = gradient_WM(x,A);
      zk = FWR_dir(x,real(s),lb,ub);
      al = 0.7/k; 
      xh=x^-0.5;
      xh2=x^0.5;
      x  = xh2*(xh*zk*xh)^al *xh2;
      
      toc_rec=toc_rec+toc;
      
      % record results
      output(k)=wmobj(x,A);
      time(k)=toc_rec;
      grad(k)=norm(s,'fro');
      
   end
end
