function [ U, A ] = init_u (A)
	 
	 ## usage: [ U, A ] = init_u (A)
	 ## 
	 ## initialize unitary matrix U using single 
         ## matrix diagonalizer and apply
         ## this transformation to a set A
	 ##

  matsz = size(A,1);
  nmat = size(A,2)/matsz;
  
  ## find a diagonalizer by congruence (A = Q L Q^T)

  AA = A(:, 1 : matsz )*A(:, matsz + 1:2*matsz)^(-1);
  [U d] = eig(AA);

  ## transform by congruence
  for j=1:nmat
      A(:, ((j-1)*matsz + 1) : j*matsz ) = U^(-1) * A(:, ((j-1)*matsz + 1) : j*matsz ) * U;
  end

endfunction
