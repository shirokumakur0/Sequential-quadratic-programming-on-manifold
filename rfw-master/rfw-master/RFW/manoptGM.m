function out=manoptGM(A,x0,mit,algo)
% manopt toolbox | modified MW

  if 0
    thisfilefull = mfilename('fullpath');
    thisfilename = mfilename();
    thisfilepath = thisfilefull(1:strfind(thisfilefull, thisfilename)-1);

    % Add that path to the Matlab path.
    addpath(thisfilepath);

    %% Import all subpackages of Manopt

    import manopt.*;
    import manopt.tools.*;
    import manopt.solvers.*;
    import manopt.manifolds.*;

    import manopt.solvers.linesearch.*;
    import manopt.solvers.neldermead.*;
    import manopt.solvers.pso.*;
    import manopt.solvers.steepestdescent.*;
    import manopt.solvers.conjugategradient.*;
    import manopt.solvers.trustregions.*;
    import manopt.solvers.mbfgs.*;
    import manopt.manifolds.complexcircle.*;
    import manopt.manifolds.psd.*;
    import manopt.manifolds.euclidean.*;
    import manopt.manifolds.fixedrank.*;
    import manopt.manifolds.grassmann.*;
    import manopt.manifolds.oblique.*;
    import manopt.manifolds.rotations.*;
    import manopt.manifolds.sphere.*;
    import manopt.manifolds.stiefel.*;
    import manopt.manifolds.symfixedrank.*;
  end
  
  switch algo
    case 1
      method='CG';
    case 2
      method='SD';
    case 3
      method='LBFGS';
    case 4
      method='TR';
    case 5
      method='BB';
  end

  global tocs fvals counter
  counter=1;
  tocs=zeros(1,mit);
  fvals=zeros(1,mit);

  if nargin<2
    disp('Not Enough Input Arguments');
    return;
  end
  options=struct;
  if ~isfield(options,'TolX')
    options.tolgradnorm=1e-10;
  end
  if ~isfield(options,'MaxIter')
    options.maxiter=mit;
  end
  if ~isfield(options,'Display')
    options.Display='iter';
  end

  m = numel(A);
  n = size(A{1},1);

  % Create the problem structure.
  manifold=sympositivedefinitefactory(n);
  problem.M = manifold;

  % Define the problem cost function and its gradient.
  problem.costgrad = @(x) GMFGX(x,A,manifold);

  warning('off', 'manopt:getHessian:approx');
  disp(['Method = ' method]);

  tic;                                 % start keeping time

  switch method
    case 'SD', 
      [xout out.cost info] = steepestdescent(problem,x0,options);
    case 'LBFGS',
      [xout out.cost info] = rlbfgs(problem,x0,options);
    case 'TR',
      [xout out.cost info] = trustregions(problem,x0,options);  
    case 'CG',
      [xout out.cost info] = conjugategradient(problem,x0,options);
    case 'BB',
      [xout out.cost info] = barzilaiborwein(problem, x0, options);
  end
  nit=numel(info);
  for i=1:nit, out.tm(i)=info(i).time; out.it(i)=info(i).cost; end
  out.GM=xout;
end

function [f,gf]=GMFGX(x,A,manifold)
  global tocs fvals counter
  if ~isposdef(A)
    f = Inf;
    gf = Inf(size(x));
  else
    f = gmobj(x,A);
    fvals(counter) = f;
    tocs(counter) = toc;
    counter = counter+1;
    if nargout>1
      g = gradient(x,A);
      gf = manifold.egrad2rgrad(x, g);
    end
  end
  %counter
end

function g = gradient(x,A)
  g = zeros(size(x));
  for i=1:numel(A)
    g = g + logm(A{i}\x);
  end

  g=2*g/x;
  g=(g+g')/2;
end
