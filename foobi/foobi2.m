function [A, O] = foobi2 (T,isym,emtresh)
	 
	 ## usage: [U, lambda] = foobi (T)
	 ## 
	 ## Decompose 4-order tensor (Dirak symmetries)
	 
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
  elseif ( isym == 1 )
    [U, D] = tfac(C);
  endif

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

  fprintf(stdout, 'foobi2: The rank is %d\n', nvec);

  H =  U(:,1:nvec) * diag(sqrt( Ds(1:nvec) ));
#  H =  U(:,1:nvec) * diag( Ds(1:nvec) );

  H = norm_herm(H);
  
  B = formB(H);

  %% pre-whitening

  %% joint off-diagonalization
  [Q D] = joint_offdiag(B,1e-8);
  
  F = H * Q;
  [A O] = decompF(F,isym);

endfunction
