function M = unpackM (R)
	 
	 ## usage: W = unpackM (R)
	 ## 
	 ## Unpacks upper triangular matrices, stored as 
	 ## columns of R to a 3-order tensor T_ijk.
	 ## The order of elements in a vector R_k of R is:
	 ## R_k = [ w_11 .. w_rr w_12 .. w_1r .. w_r-1,r ]
	 ## The order of entries in T is  
	 ## 
	 ##
	 ## As a remark, the size of R is r(r+1)/2 x r

  nvec   = size(R,2);
  
  M  = zeros(nvec, nvec, nvec);
  mj = zeros(nvec);
  
  for j = 1:nvec

    %% unpack diagonal 
    mj = diag(R(1:nvec,j));
    pos = nvec + 1;

    %% unpack upper triangle
    for k = 1:nvec-1
      mj(k,k + 1:nvec) = R(pos:pos + nvec - 1 - k,j);
      pos += nvec - k;
    end

    %% restore full symmetric matrix
    mj = 1/2 * (mj + mj.');

    %% stack
    M(:,:,j) = mj;
  end

endfunction
