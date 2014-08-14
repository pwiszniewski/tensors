function [Q D] = tfac (A,nSteps,verb)
	 
	 ## usage: [Q D] = tfac (A, nSteps, verb)
	 ## 
	 ## Decompose a complex symmetric matrix using Takagi 
	 ## factorization A = Q D Q^T
	 ## where D is diagonal, Q is unitary
	 ##
	 ## nSteps - maximal number of steps in the Lanczos algorithm
         ## verb   - verbosity

  % get size of matrix
  if (~ issquare(A) || ~issymmetric(A,1e-8))
    error ('A is not symmetric or square');
  else
    n = size(A,1);
  endif
 
  if (~exist('nSteps','var'))
    nSteps = 100; %100 steps of Lanczos iterations by default
  endif
 
  if (~exist('verb','var'))
     verb = 0;
  endif
 
  % Lanczos tridiagonalization using modified partial orthogonalization
  % and restart
  [a,b,Q1,nSteps,nVec] = LanMPOR(A,rand(n,1),nSteps);
  
  % calculate and report errors
  if (verb > 0)
    tmp = norm(Q1'*Q1 - eye(nSteps), 'fro')/(nSteps*nSteps);
    fprintf('\nError in orthogonality in tridiagonalization: %E', tmp);
    tmp = norm(Q1'*A*conj(Q1) - (diag(a)+diag(b,1)+diag(b,-1)), 'fro');
    tmp = tmp/(nSteps*nSteps);
    fprintf('\nError in tridiagonalization: %E', tmp);
  endif

  % prompt the user as to which SVD algorithm to use
  svd_selection = 1;

  % svd of tridiagonal 
  if (svd_selection == 1)
    [s,Q2] = CSSVD(a, b);		% pure QR 
  elseif (svd_selection == 2)
    [s,ifail,Q2] = cstsvdd(a, b);	% divide-and-conquer
  else
    [s,Q2] = cstsvdt(a, b);             % twisted factorization
  end

  Q = Q1*Q2;
  D = diag(s);
  % check results
  if (verb > 0)
    tmp = norm(Q'*Q - eye(nSteps), 'fro')/(nSteps*nSteps);
    fprintf('\nError in orthogonality: %E', tmp);

    % calculate and report errors
    if nSteps==n

      tmp = norm(A - Q*diag(s)*conj(Q'), 'fro')/(n*n);
      fprintf('\nError in Takagi factorization: %E\n', tmp);
    end
  endif
	 
endfunction
