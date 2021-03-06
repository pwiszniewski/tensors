function [ U D ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 
	 ## usage: [ U A ] = jacobi (A, delta, maxcyc, optgoal, verb)
	 ## 
	 ## jacobi iteration for a set of matrices
	 ##

  if (ndims(A) != 2)
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
  
  [indx, B] = init_indices(A);
  [normf, other_normf] = get_norms(B, optgoal);
  
  %check for early exit
  if (other_normf < 1e-6)
    U = eye(matsz);
    D = A;
    fprintf(stdout, 'Diagonal matrix. jacobi exited early\n');    
    return
  end
     
  ## pre-whitening
  [U, A] = init_u(A);

  if (verb > 0) 
    fprintf(stdout, 'Starting Jacobi iteration. Initial norm %19.4e, off-norm: %19.4e\n', normf, other_normf);
  endif  

  other_normf_old = other_normf + 3*delta;

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
    pair = [k;l];
    indk = k : matsz : matsz*nmat;
    indl = l : matsz : matsz*nmat;

    h = form_h(A,k,l);
    
    [c,s] = find_angle(h,optgoal);
    
    R = [c -conj(s) ; s c ];
    
    ## apply U * R       
    U(:,pair) = U(:,pair) * R;
    
    ## apply R' * A * R
    A(pair,:) = R' * A(pair,:);
    A(:,[indk indl])    = [ c*A(:,indk)+s*A(:,indl) \
			-conj(s)*A(:,indk)+c*A(:,indl) ] ;

    [indx, B] = update_indices(indx, B, A, k, l);
    [normff, other_normff] = get_norms(B, optgoal);
    
    if (verb > 1) 
       fprintf(stdout, 'sweep: %d, current off-norm: %19.4e\n', it, other_normf);
    endif
  end

  if (verb > 0) 
    fprintf(stdout, 'Finished Jacobi iteration. Sweeps: %d, final off-norm: %19.4e\n', it, other_normf);
  endif  

  D = A;

endfunction
