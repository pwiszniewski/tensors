function  normof = bfgs_wrapper (cvec, dim, EI, isym, normuse)
	 
	 ## usage: normof = bfgs_wrapper (cvec, dim, EI)
	 ## 
	 ## 
	  
  if (length(cvec) == dim*(dim-1)/2)
    S = vec2triu(cvec,dim,0);
    d = diag(S);
    S = S + S.';
    S = S - d;
  elseif (length(cvec) == dim*dim)
    S = vec2triu( cvec(1 : dim*(dim+1)/2),dim,0);
    d = diag(S);
    S = S + S';
    S = S - diag(d);
    A = vec2triu(cvec(dim*(dim+1)/2 + 1:dim*dim),dim,1);
    A = A - A.';
  else
    error('wrong length of the parameter vector %d\n', length(cvec));
  end
	 
  U = expm(i*S - A);

  EI_new = t2e(EI,U,isym);
  
  if ( normuse == 1 )
    [ normd normof ] = diagnorm(EI_new,isym);
  else
    [ normd normof ] = bidiagnorm(EI_new,isym);
  endif

endfunction
