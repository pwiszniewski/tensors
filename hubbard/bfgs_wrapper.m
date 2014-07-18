function  normof = bfgs_wrapper (cvec, dim, EI, normuse)
	 
	 ## usage: normof = bfgs_wrapper (cvec, dim, EI)
	 ## 
	 ## 
  packing = 0; %mulliken

  H = reshape(cvec,dim,dim);
  S = triu(H) + triu(H,1)';
  A = - tril(H,-1) + tril(H,-1)';
  U = expm(i*S - A);

  EI_new = t2e(EI,U,packing);
  
  if ( normuse == 1 )
    [ normd normof ] = diagnorm(EI_new);
  else
    [ normd normof ] = bidiagnorm(EI_new);
  endif

endfunction
