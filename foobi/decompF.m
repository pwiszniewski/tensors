function [A O] = decompF (F)
	 
	 ## usage: [A O] = decompF (F)
	 ## 
	 ## This function is a final step of the FOOBI
	 ## algorithms. The columns a_k of the factor matrix A
	 ## are acessed as the left singular vectors of the 
	 ## matrices f_k , which are stacked as columns of F.
	 ##
	 ## The amplitudes O are the singular values associated 
	 ## with A_j (this is a hypotesis and is not explicitly 
	 ## tested).

  veclen = size(F,1);
  matsz  = sqrt(veclen);
  nvec   = size(F,2);
  
  A = zeros(matsz,nvec);
  O = zeros(nvec,1);

  for j = 1:nvec
     [u s v]  = svds ( reshape(F(:,j),matsz,matsz), 1, 'L' );
     A(:,j) = u;
     O(j) = s;
  end

endfunction
