function A = decompF (F)
	 
	 ## usage: A = decompF (F)
	 ## 
	 ## This function is a final step of the FOOBI
	 ## algorithms. The columns a_k of the factor matrix A
	 ## are acessed as the left singular vectors of the 
	 ## matrices f_k , which are stacked as columns of F.
	 
  veclen = size(F,1);
  matsz  = sqrt(veclen);
  nvec   = size(F,2);
  
  A = zeros(matsz,nvec);

  for j = 1:nvec
     [u s v]  = svds ( reshape(F(:,j),matsz,matsz), 1, 'L' );
     A(:,j) = u;
  end

endfunction
