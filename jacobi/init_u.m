function [ U, A ] = init_u (A)
	 
	 ## usage: [ U, A ] = init_u (A)
	 ## 
	 ## initialize unitary matrix U using single 
         ## matrix diagonalizer and apply
         ## this transformation to a set A
	 ##

	 
  a_herm = 1/2 * (A(:,:,1) + A(:,:,1)');

  ## returns right eigenvectors
  [U d] = eig(a_herm);
  
  U = U';

  nmat = size(A,3);
  
  for j=1:nmat
      A(:,:,j) = U * A(:,:,j) * U';
  end

endfunction
