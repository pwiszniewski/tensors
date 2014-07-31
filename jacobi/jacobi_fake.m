function [ Q D ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 
	 ## usage: [ Q A ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 ## 
	 ## jacobi iteration for a set of matrices
	 ##

  if (ndims(A) != 2)
     error('Data is of wrongly shaped!');
  end

  matsz = size(A,1);
  nmat = size(A,2)/matsz;

  if (~exist('delta','var'))
     delta = 1e-8;
  endif

  if (~exist('maxcyc','var'))
     maxcyc = 100;
  endif

  if (~exist('optgoal','var'))
     optgoal = 0;   ## diagonalization. zero-diagonalization logic not impl
  endif

  if (~exist('verb','var'))
     verb = 1;  
  endif

  ## pre-whitening
  
  AA = A(:,1:matsz) * A(:,matsz+1:2*matsz)^(-1);

  ## left and right singular vectors
  [Q d] = eig(AA);

  for j=1:nmat
      A(:, ((j-1)*matsz + 1) : j*matsz ) = Q^(-1) * A(:, ((j-1)*matsz + 1) : j*matsz ) * (Q.')^(-1);
  end

  D = A;

endfunction
