function [A O] = foobi (T,isym,emtresh)
	 
	 ## usage: [A, O] = foobi (T, isym, emtresh)
	 ## 
	 ## Decompose 4-order tensor (Dirak symmetries)
	 ## isym - symmetry of the input tensor
	 ##      = 0 (default) - cumulant
         ##      = 1 Mulliken
	 ##      = 2 Dirak

  dimin = size(T,4);

  if (~exist('isym','var'))
     isym = 0;
     warning('using cumulant symmetry by default');
  endif 

  if (~exist('emtresh','var'))
    emtresh = 1e-8;
  endif

  C = reshape(T,dimin*dimin, dimin*dimin);

  if ( (isym == 0) || (isym == 2) )
    [U, D, ~] = svd(C);
##    [U, D] = eig(C);
  elseif ( isym == 1 )
    [U, D] = tfac(C);
  endif

  %% Warning! have to check that eigenvalues are real each time
  [Ds indD] = sort(real(diag(D)),'descend'); 
  U = U(:,indD);
  
  %% truncation of non-significant eigenvalues and eigenvectors
  nvec = 1;
  for j = 1:dimin*dimin
      if ( abs(Ds(j)) < emtresh ) 
	 nvec = j - 1;
	 break
      endif
      nvec = j;
  end

  fprintf(stdout, 'foobi: The rank is %d\n', nvec);
  fprintf(stdout, 'foobi: the lowest kept eigenvalue: %f\n', Ds(nvec));
  fflush(stdout);

  H =  U(:,1:nvec) * diag(sqrt( Ds(1:nvec) ));
#  H =  U(:,1:nvec) * diag( Ds(1:nvec) );

  H = norm_herm(H);

#  P = formP_full(H); %% check this!
  P = formP(H,isym);
  
  %% since we do full svd we take only nvec right sing. vectors

  [~, ~, R] = svd(P);
  R = R(:,end:-1:end - nvec + 1);

  #W = unpacktri(R(:,1:nvec));

  M  = unpackM(R);
  [Q out] = cpd3_sgsd(M,{eye(nvec),eye(nvec),eye(nvec)},struct('TolFun',1e-6,'MaxIter',5000));

  if ( abs(out.fval(end)) > 1e-6 )
     warning ('foobi: convergence problem in joint diagonalization, conv = %e, niter = %d\n', out.fval(end), length(out.fval));
  end

  F = H * Q{1};
  [A O] = decompF(F,isym);

endfunction
