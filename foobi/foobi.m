function [A O D] = foobi (T,isym,nvemax,emtresh)
	 
	 ## usage: [U, lambda] = foobi (T)
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

  if (~exist('nvemax','var'))
    nvemax = 6;
  endif
  
  if (nvemax > dimin*dimin)
     nvemax = dimin*dimin;
     warning('setting nvemax to maximal possible rank %d', dimin*dimin);
  endif 
  
  C = reshape(T,dimin*dimin, dimin*dimin);

  if ( (isym == 0) || (isym == 2) )
    [U D] = eig(C);
  else
    [U D UU] = svd(C);
  end

  [Ds indD] = sort(real(diag(D)),'descend'); 
  U = U(:,indD);
  
  %% truncation of non-significant eigenvalues and eigenvectors
  nvec = 1;
  for j = 1:nvemax
      if ( abs(Ds(j)) < emtresh ) 
	 nvec = j - 1;
	 break
      endif
      nvec = j;
  end

  %% H =  U(:,1:nvec) * diag(sqrt( Ds(1:nvec) ));
  H =  U(:,1:nvec) * diag( Ds(1:nvec) );

  H = norm_herm(H);

  P = formP(H); %% check this! this is wrong for cumulants

  %% the following may cause problems if the eigenvalues are quite small
  %% need an idea on how to choose an appropriate guess.
  %[L S R] = svds(P,nvec,1e-6); 

  %% since we do full svd we take only nvec right sing. vectors

  [L S R] = svd(P);

  [Ss indS] = sort(diag(S),'ascend'); 
  R = R(:,indS);

  W = unpacktri(R(:,1:nvec));

#  [Q Di] = jacobi_fake(W,1e-8);

  [Q L] = jacobi(W,1e-8);
#  [Q L] = joint_diag(W,1e-8);

  F = H * Q;
  [A O] = decompF(F);

endfunction
