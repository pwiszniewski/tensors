function W = unpacktri (R)
	 
	 ## usage: W = unpacktri (R)
	 ## 
	 ## Unpacks upper triangular matrices, stored as 
	 ## columns of R to square, which are stored stacked in W.
	 ## The order of elements in a vector R_k of R is:
	 ## R_k = [ w_11 .. w_rr w_12 .. w_1r .. w_r-1,r ]
	 ## The order of matrices W_k in W is 
	 ## W = [ W_1 W_2 .. W_R ]
	 ##
	 ## As a remark, the size of R is r(r+1)/2 x r

  nvec   = size(R,2);
  
  W  = zeros(nvec, nvec*nvec);
  wj = zeros(nvec);
  
  for j = 1:nvec

    %% unpack diagonal 
    wj = diag(R(1:nvec,j));
    pos = nvec + 1;

    %% unpack upper triangle
    for k = 1:nvec-1
      wj(k,k + 1:nvec) = R(pos:pos + nvec - 1 - k,j);
      pos += nvec - k;
    end

    %% restore full symmetric matrix
    wj = 1/2 * (wj + wj.');

    %% stack
    ind = [1+(j-1)*nvec : j*nvec];
    W(:,ind) = wj;
  end

endfunction
