function [ U, A ] = init_u (A)
	 
	 ## usage: [ U, A ] = init_u (A)
	 ## 
	 ## initialize unitary matrix U using single 
         ## matrix diagonalizer and apply
         ## this transformation to a set A
	 ##

  matsz = size(A,1);
  nmat = size(A,2)/matsz;

  a_herm = 1/2 * (A(:,1:matsz) + A(:,1:matsz)');

  ## returns right eigenvectors
  [U d] = eig(a_herm);
  
  for j=1:nmat
      A(:, ((j-1)*matsz + 1) : j*matsz ) = U' * A(:, ((j-1)*matsz + 1) : j*matsz ) * U;
  end

endfunction
