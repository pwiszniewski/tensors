function [ EI, U Z ] = initrandgen (dim, u_cplx, z_cplx, z_mode, packing)
	 
	 ## usage: EI U Z = initrandrec (dim, u_cplx, z_cplx, z_mode, packing)
	 ## Builds integrals from rectangular factors, which was used
         ## in older numerical experiments
	 ## 
	 ##  u_cplx = 1 - use complex unitary factors 
	 ##           0 - use real factors
	 ##  z_cplx = 1 - use coplex Z
	 ##           0 - use real Z
         ##
	 ##  z_mode = 1 - identity
	 ##           2 - diagonal
	 ##           3 - symmetric
	 ##           4 - antisymmetric
	 ##           5 - hermitian     (complex only)
	 ##           6 - antihermitian (complex only)
	 ##           7 - orthogonal/unitary
	 ## packing = 1 Dirak
	 ##           0 Mulliken 

	 
  if (z_cplx != 1)
    if (z_mode == 1)
      Z = eye(dim*(dim+1)/2);  ## Real evals, but not ones | real
      ## evals, not ones
    elseif (z_mode == 2)
      Z = diag(rand(dim*(dim+1)/2,1));	## Real evals | real evals   
    elseif (z_mode == 3)
      Z = rand(dim*(dim+1)/2); ## complex paired evals, n(n+1)/2 total
      ## | real evals
      Z = 1/2 * (Z + Z.');
    elseif (z_mode == 4)
      Z = rand(dim*(dim+1)/2); ## breaks symm, 5 | breaks symm. 6
      Z = 1/2 * (Z - Z.');
    elseif (z_mode == 7)
      Z = orth(rand(dim*(dim+1)/2)); ## breaks symm, 5 | breaks symm. 6
    elseif (z_mode == 8)
      Z = rand(dim*(dim+1)/2); ## breaks symm, 5 | breaks symm. 6
    else
      error('initrandrec', 'z_mode not supported\n');
    end
  else
    if (z_mode == 1)
      Z = eye(dim*(dim+1)/2) + i*eye(dim*(dim+1)/2); # breaks symm, 6
				# | breaks symm. 6
    elseif (z_mode == 2)
      Z = diag(rand(dim*(dim+1)/2,1) + i*rand(dim*(dim+1)/2,1));    ##
      ##breaks symm, 6 | breaks symm. 6
    elseif (z_mode == 3)
      Z = rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2); ## breaks symm,
      ## 6 (nondiagonalizable) | breaks symm 6
      Z = 1/2 * (Z + Z.');
    elseif (z_mode == 4)
      Z = rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2); ## breaks symm,
      ## 7 | breaks symm 6
      Z = 1/2 * (Z - Z.');
    elseif (z_mode == 5)
      Z = rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2); ## breaks symm,
      ## 3 | real evals
      Z = 1/2 * (Z + Z');
    elseif (z_mode == 6)
      Z = rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2); ## breaks symm,
      ## 7 | breaks symm 6
      Z = 1/2 * (Z - Z');
    elseif (z_mode == 7)
      Z = orth(rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2)); ## breaks
      ## symm, 7 | breaks symm 6 
    elseif (z_mode == 8)
      Z = rand(dim*(dim+1)/2) + i*rand(dim*(dim+1)/2); ## breaks symm,
      ## 7 | breaks symm 6
    else
      error('initrandrec', 'z_mode not supported\n');
    end
  end

  if (u_cplx != 1)
    U = rand(dim*(dim+1)/2, dim);
    U = orth(U);
  else
    U = rand(dim*(dim+1)/2, dim) + i*rand(dim*(dim+1)/2, dim);
    U = orth(U);
  endif

 EI = zeros(dim,dim,dim,dim);

 if (packing == 0)
   for j = 1 : dim
     for k = 1 : dim
       for l = 1 : dim
	 for m = 1 : dim
	   for p = 1 : (dim*(dim+1)/2)
	     for q = 1 : (dim*(dim+1)/2)
	       EI(j,k,l,m) += conj(U(p,j)) * U(p,k) * Z(p,q) * conj(U(q,l)) \
			      * U(q,m);
	     end
	   end
	 end
       end
     end
   end
   
 else
   for j = 1 : dim
     for k = 1 : dim
       for l = 1 : dim
	 for m = 1 : dim
	   for p = 1 : (dim*(dim+1)/2)
	     for q = 1 : (dim*(dim+1)/2)
	       EI(j,k,l,m) += conj(U(p,j)) * conj(U(p,k)) * Z(p,q) * U(q,l) \
			      * U(q,m);
	     end
	   end
	 end
       end
     end
   end

 end
     
endfunction
