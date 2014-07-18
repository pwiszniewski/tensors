function [ U A ] = jacobi_off (A, delta, maxcyc, optgoal, verb)
	 
	 ## usage: [ U A ] = jacobi_off (A, delta, maxcyc, optgoal, verb)
	 ## 
	 ## jacobi iteration for a set of matrices
	 ##

  nmat = size(A,3);
  matsz = size(A,1);

  if (~exist('delta','var'))
     delta = 1e-8;
  endif

  if (~exist('maxcyc','var'))
     maxcyc = 100;
  endif

  if (~exist('optgoal','var'))
     optgoal = 1;
  endif

  if (~exist('verb','var'))
     verb = 1;  
  endif

  ## pre-whitening step
  #[U, A] = init_u(A);
  U = eye(matsz);
  
  [indx, B] = init_indices(A, optgoal);
  [normf, other_normf] = get_norms(B, optgoal);
  other_normf_old = other_normf + 3*delta;

  if (verb > 0) 
    fprintf(stdout, 'Starting Jacobi iteration. initial normf: %19.4e\n', other_normf);
  endif
  
  
  it = 0;
  while (  (it < maxcyc))
#(abs(other_normf - other_normf_old) > delta) &&
    it += 1;

#    k = 1;
#    for j = 2:matsz-1
#      if ( B(j,indx(j)) < B(k,indx(k)))
#	 k = j;
#      endif
#    end
#    l = indx(k);

  for k = 1:matsz-1
    for l = k + 1:matsz

    h = form_h(A,k,l);
    
    [c,s] = find_angle(h,optgoal);
    
    ## apply R * A * R^h
    A = rotate_set(A,k,l,c,s);
    
    ## apply R * U       
    for j = 1:matsz
      [U(k,j), U(l,j)] = rotate(U(k,j),U(l,j),c,s);
    end
    
    [indx, B, diffB] = update_indices(indx, B, A, k, l, optgoal);
    # [normff, other_normff] = get_norms(B, optgoal);
    
    other_normf_old = other_normf;
    other_normf = update_norm(other_normf, diffB, k, l, optgoal);

    end
  end

    if (verb > 0) 
       fprintf(stdout, 'sweep: %d, residual norm: %19.4e\n', it, other_normf);
    endif
  end

endfunction
