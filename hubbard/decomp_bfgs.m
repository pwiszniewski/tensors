function [ U_init U EI ] = decomp_bfgs (dim, cplx, guess, normuse)
	 
	 ## usage: [ normd normof ] = decomp_bfgs (dim, cplx, guess, normuse)
	 ## 
	 ## guess   = 1 - read
	 #         <> 1 - new 
	 ## normuse = 1 - diagnorm
	 ##        <> 1 - bidiagnorm

  isym = 1; % mulliken, for dirak change bidiagnorm

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
    [ EI U_init ] = initrandhub(dim, cplx, isym);  

    if (cplx == 1)
      cvec = rand(dim*dim,1);;
    else
      cvec = rand(dim*(dim-1)/2);
    endif

  endif

  args = {cvec, dim, EI, isym, normuse};

  [ cvecmin, minnorm, conv, iters ] = bfgsmin("bfgs_wrapper", args,  control );

#  fprintf("Results dumped to optresults\n");
#  save optresults EI cvecmin;

  if (cplx == 1)
    S = vec2triu(cvec(1:dim*(dim+1)/2),dim,0);
    d = diag(S);
    S = S + S';
    S = S - diag(d);
    A = vec2triu(cvec(dim*(dim+1)/2 + 1:dim*dim),dim,1);
    A = A - A.';
  else
    S = vec2triu(cvec,dim,0);
    d = diag(S);
    S = S + S.';
    S = S - d;
  endif    

  U = expm(i*S - A);

  EI = t2e(EI,U,isym);
  
endfunction
