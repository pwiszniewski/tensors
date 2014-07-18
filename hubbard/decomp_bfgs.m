function decomp_bfgs (dim, cplx, guess, normuse)
	 
	 ## usage: [ normd normof ] = decomp_bfgs (dim, cplx, guess, normuse)
	 ## 
	 ## guess   = 1 - read
	 #         <> 1 - new 
	 ## normuse = 1 - diagnorm
	 ##        <> 1 - bidiagnorm

  packing = 0; % mulliken, for dirak change bidiagnorm

  maxit = 10000;
  verb = 2;
  convcrit = 1;
  tol = 1e-6;
  control = {maxit, verb, convcrit};

  if ( guess == 1 )
    A = load('optresults');
    EI = A.EI;
    cvec = A.cvecmin;
  else
    EI = initrandei(dim, cplx, packing);  
    cvec = vec(rand(dim));
  endif

  args = {cvec, dim, EI, packing, normuse};

  [ cvecmin, minnorm, conv, iters ] = bfgsmin("bfgs_wrapper", args,  control );

  fprintf("Results dumped to optresults\n");
  save optresults EI cvecmin;

endfunction
