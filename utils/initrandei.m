function EI = initrandei (dimin, rnk, isym, cplx)
	 
	 ## usage: EI = initrandei (dimin, rnk, isym, cplx)
	 ## returns random integrals in the Mulliken or Dirak ordering 
	 ##  
	 ## isym = 1 Mulliken
         ##        2 Dirak
	 ##      

  if ( cplx == 1 )
    EI = rand(dimin,dimin,dimin,dimin) + i*rand(dimin,dimin,dimin,dimin);
  else
    EI = rand(dimin,dimin,dimin,dimin);
  endif

  EI = symmetrize(EI,isym,cplx);

  ET = reshape(EI,dimin*dimin, dimin*dimin);
  [V S VV] = svd(ET);
  ss = diag(S);
  ss(end:-1:rnk+1) = 0;
  ET = V * diag(ss) * VV';
  EI = reshape(ET,dimin,dimin,dimin,dimin);

endfunction
