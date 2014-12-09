function [EI trunc ] = trunctensvd (ET,rnk)

	 ## usage: [EI trunc]= trunctensvd (EI, rnk)
	 ## returns alternated 4-tensor by truncating an SVD of it's 
         ## unfolding
	 ##  

  dim = size(ET,1);
  [U S V] = svd(reshape(ET,dim*dim,dim*dim));
  ss = diag(S);
  trunc = ss(rnk+1);
  ss(end:-1:rnk+1) = 0;
  ET = U * diag(ss) * V';
  EI = reshape(ET,dim,dim,dim,dim);
  
endfunction
