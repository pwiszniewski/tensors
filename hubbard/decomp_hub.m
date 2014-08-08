function [ EI U_init U ] = decomp_hub (dim, cplx)
	 
	 ## usage: [ EI_new U_init U ] = decomp_hub (dim, cplx)
	 ## 
	 ## Decompose random 4-index tensor of dimension dim x dim x
         ## dim x dim and rank dim
	 ## into unitary matrices
	 ##
	 
  A = rand(dim);
  B = rand(dim);

  if ( cplx == 1 )
     A = A + i*rand(dim);
     B = B + i*rand(dim);
  endif

  FA = zeros(dim);
  FB = zeros(dim);

  [ EI U_init ] = initrandhub(dim, cplx, 1);

  ## Warning! This contraction is only valid for packing = Mulliken

  for j = 1 : dim
    for k = 1 : dim
	FA(j,k) = trace(vec(EI(j,k,:,:))*transpose(vec(A)));
	FB(j,k) = trace(vec(EI(j,k,:,:))*transpose(vec(B)));	
    end
  end

 [ U D ] = eig(FA,FB); 
 S = U'*U;
 
 ## normalize 
 for j = 1 : dim
   U(:, j) = U(:, j)./ sqrt(S(j,j));
 end

 ## transform using U

 EI_new = t2e(EI, U, 1);

 ## offdiagonal norm should be zero, diagonal norm should be 8 x dim 
 
 [ normd_new normof_new ] = bidiagnorm(EI_new, 1)

 ## should be complex matrix with unit columns
 
 U_init * U;

endfunction


