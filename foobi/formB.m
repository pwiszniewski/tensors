function B = formB (H)
	 
	 ## usage: B = formB (H)
	 ## 
	 ## Form the matrices B_ij, which should be off 
	 ## diagonalized in the FOOBI2 method. 
	 ## B is stacked column wize as follows
	 ## B = [ real(b_1) .. real(b_n) imag(b_1) .. imag(b_n) ]
	 ##
	 ## B_ij is defined as ( B_ij )_st = ( phi(H_s, H_t) )_ij
	 ##
	 ## phi(X,Y) = XY - trace(X)Y + YX - trace(Y)X 

  veclen = size(H,1);
  matsz  = sqrt(veclen);
  nvec   = size(H,2);

  B = zeros(nvec , nvec * matsz * matsz * 2);
  # temporary array

  for s = 1:nvec
    Hs = reshape(H(:,s),matsz,matsz);

    for t = 1:nvec
	Ht = reshape(H(:,t),matsz,matsz);
	PHI = Hs * Ht + Ht * Hs - trace(Hs) * Ht - \
	      trace(Ht) * Hs;
	indr = [t : nvec : nvec * matsz * matsz];
	indi = [nvec * matsz * matsz + t : nvec : nvec * matsz * \
						  matsz * 2];
	B(s,indr) = real(vec(PHI));
	B(s,indi) = imag(vec(PHI));
    end
  end
	 
endfunction
