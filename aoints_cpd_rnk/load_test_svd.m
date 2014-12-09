function [nbf, nval, res] = load_test_svd (cutoff)
	 
	 ## usage: [nbf, nval, res] = load_test_svd (cutoff)
	 ## 
	 ## Script to load integrals from "aoints/integrals.h5",
	 ## do an SVD of them, cut anything lower than cutoff
	 ##
	 ## nbf  - number of basis functions
	 ## nval - number of large SVD values
	 ## res  - smallest value kept 
	 ##
	 
	 h = load("aoints/integrals.h5");
	 nbf = h.integrals.Nbf(1);
	 [~, s, ~]  = svd(sparse(reshape(h.integrals.TwoBody,nbf*nbf, nbf*nbf)));
	 
	 s = diag(s);
	 
	 nval = 1;
	 for j = nbf*nbf:-1:1
	   if ( abs(s(j)) > cutoff ) 
	     nval = j;
	     break
	   endif
	 end

	 res = s(nval);
	 
endfunction
