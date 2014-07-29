function [ U D ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 
	 ## usage: [ U A ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 ## 
	 ## jacobi iteration for a set of matrices
	 ##

  if (ndims(b) != 2)
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

  ## pre-whitening step
  [U, A] = init_u(A);
  
  [indx, B] = init_indices(A)
  [normf, other_normf] = get_norms(B, optgoal);
  other_normf_old = other_normf + 3*delta;

  if (verb > 0) 
    fprintf(stdout, 'Starting Jacobi iteration. Initial norm %19.4e, off-norm: %19.4e\n', normf, other_normf);
  endif  
  
  it = 0;
  while ( (abs(other_normf - other_normf_old) > delta) && (it < maxcyc))
    it += 1;

    k = 1;
    for j = 2:matsz-1
      if ( B(j,indx(j)) > B(k,indx(k)))
	 k = j;
      endif
    end
    l = indx(k);

    h = form_h(A,k,l);
    
    [c,s] = find_angle(h,optgoal);
    
    ## apply R * A * R^h
    A = rotate_set(A,k,l,c,s);
    
    ## apply R * U       
    for j = 1:matsz
      [U(k,j), U(l,j)] = rotate(U(k,j),U(l,j),c,s);
    end
    
    [indx, B, diffB] = update_indices(indx, B, A, k, l);
    # [normff, other_normff] = get_norms(B, optgoal);
    
    other_normf_old = other_normf;
    other_normf = update_norm(other_normf, diffB, k, l, optgoal);
    
    if (verb > 1) 
       fprintf(stdout, 'sweep: %d, current off-norm: %19.4e\n', it, other_normf);
    endif
  end

  if (verb > 0) 
    fprintf(stdout, 'Finished Jacobi iteration. Sweeps: %d, final off-norm: %19.4e\n', it, other_normf);
  endif  

  D = reshape(A,matsz,matsz*nmat);

endfunction
