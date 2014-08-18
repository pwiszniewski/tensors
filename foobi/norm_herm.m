function H = norm_herm (H)
	 
	 ## usage: H = norm_herm (H)
	 ## 
	 ## Hermitize or symmetrize the eigenmatrices h_k stored 
         ## as column vectors in H. Initially, h_k may not be hermitean
         ## or complex symmetric. 
         ## To rotate it to the hermitean (or symmetric) matrix 
	 ## we need to multiply it by alpha*I. alpha is built as
         ## follows:
	 ## if any diagonal element Z of h_k is nonzero, 
	 ## alpha is a unit complex factor that makes Z real.
	 ## if all diagonal elements are zero, than we should 
	 ## choose a phase such that (h_k)_ij = (h_k)_ji (not impl).
	 
  veclen = size(H,1);
  matsz  = sqrt(veclen);
  nvecs  = size(H,2);
  norm_vec = zeros(nvecs,1);
  
  for j = 1:nvecs
    for k = 1:matsz + 1:veclen
	if ( abs( H(k,j) ) != 0 )
	  norm_vec(j) = conj( H(k,j) / abs(H(k,j)) );

	  %% avoid multiplication by negative numbers, as it changes the
	  %% direction of eigenvectors (POSSIBLY UNNECESSARY)
	  if ( real(norm_vec(j)) < 0 )
	     norm_vec(j) = -norm_vec(j);
	  endif

	  break

	endif
    end
  end

  H = H * diag(norm_vec);

endfunction
