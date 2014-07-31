function P = formP_test (H)
	 
	 ## usage: P = formP (H)
	 ## 
	 ## Forms a tensors P_j used in the FOOBI 
         ## algorithm. Each column of P is a vectorized
	 ## 4-index tensor, built of the matrices h_k, which 
	 ## are stored vectorized as columns of H. 
         ## If H_s, H_t are some two matrices stored in H, than  
         ## a column P_st is defined as  
	 ## (P_st)_ijkl = (H_s)_ij * (H_t)_kl + (H_t)_ij * (H_s)_kl -
	 ## -(H_s)_ik * (H_t)_jl - (H_t)_ik * (H_t)_jl
	 ## 
	 ## The order of entries in P is [ P_11 P_22 .. P_RR P_12 ..
         ## P_1R .. P_R-1R ] (only diagonal and upper triangle is
         ## stored, as is allowed by symmetry of P)
	 
  veclen = size(H,1);
  dimin  = sqrt(veclen);
  nvec   = size(H,2);
  
  P = zeros(nvec,nvec,veclen);

  for j = 1:nvec
    hj = reshape(H(:,j),dimin,dimin);      
    P(:,j) = 2 * ( vec(kron(hj.',hj)) - vec(H(:,j)*H(:,j)') );
  end

  jk = nvec;
  for j = 1:nvec - 1
    hj = reshape(H(:,j),dimin,dimin);

    for k = j + 1:nvec
      jk += 1;
      hk = reshape(H(:,k),dimin,dimin);

      P(:,jk) = 2 * ( vec(kron(hk.',hj)) + vec(kron(hj.',hk)) - \
		vec(H(:,j)*H(:,k)') - vec(H(:,k)*H(:,j)') );
    end
  end
     
endfunction
